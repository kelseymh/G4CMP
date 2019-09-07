/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
//
// 20140324  Drop hard-coded IV scattering parameters; get from lattice
// 20140418  Drop valley transforms; get from lattice
// 20170802  Remove MFP calculation; use scattering-rate model
// 20190704  Add selection of rate model by name, and material specific
// 20190906  For rate model selection, pass string by value
// 20190906  Push selected rate model back to G4CMPTimeStepper for consistency

#ifndef G4CMPInterValleyScattering_h
#define G4CMPInterValleyScattering_h 1

#include "globals.hh"
#include "G4CMPVDriftProcess.hh"


class G4CMPInterValleyScattering : public G4CMPVDriftProcess { 
public:
  G4CMPInterValleyScattering();
  virtual ~G4CMPInterValleyScattering();

  // Select different rate models by string (globally or by material)
  using G4CMPVProcess::UseRateModel;		// Avoid function hiding
  void UseRateModel(G4String model);		// Non-virtual to use in ctor

  // Do scattering action here
  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

  // Only electrons have physical valleys associated with them
  virtual bool IsApplicable(const G4ParticleDefinition&);

protected:
  // Change registered scattering rate based on material, if necessary
  virtual G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*);

private:
  G4String modelName;		// Last chosen rate model, to avoid memory churn

  void PushModelToTimeStepper();	// Ensure model is used for stepping

private:
  //hide assignment operator as private
  G4CMPInterValleyScattering(G4CMPInterValleyScattering&);
  G4CMPInterValleyScattering& operator=(const G4CMPInterValleyScattering& right);
};

#endif

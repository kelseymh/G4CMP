/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef G4CMPQPScattering_h
#define G4CMPQPScattering_h 1

#include "globals.hh"
#include "G4CMPVProcess.hh"
#include "G4CMPKaplanUtils.hh"


class G4CMPQPScattering : public G4CMPVProcess, G4CMPKaplanUtils {
public:
  G4CMPQPScattering();
  virtual ~G4CMPQPScattering();

  //Using Rate models
  using G4CMPVProcess::UseRateModel;            // Avoid function hiding
  void UseRateModel(G4double model);            // Non-virtual to use in ctor

  virtual bool IsApplicable(const G4ParticleDefinition&);

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&) override;

  G4CMPQPScattering(G4CMPQPScattering&) = delete;
  G4CMPQPScattering(G4CMPQPScattering&&) = delete;
  G4CMPQPScattering& operator=(const G4CMPQPScattering&) = delete;
  G4CMPQPScattering& operator=(const G4CMPQPScattering&&) = delete;

protected:
  virtual G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*)
                                   override;


private:


};

#endif

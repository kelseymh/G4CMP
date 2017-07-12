/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
//
// 20170711  New class implementing physically motivated IV scattering

#ifndef G4CMPIVScatteringPhysical_h
#define G4CMPIVScatteringPhysical_h 1

#include "globals.hh"
#include "G4CMPVDriftProcess.hh"


class G4CMPIVScatteringPhysical : public G4CMPVDriftProcess { 
public:
  G4CMPIVScatteringPhysical();
  virtual ~G4CMPIVScatteringPhysical();

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

  // Only electrons have physical valleys associated with them
  virtual bool IsApplicable(const G4ParticleDefinition&);

protected:
  virtual G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*);

private:
  //hide assignment operator as private
  G4CMPIVScatteringPhysical(G4CMPIVScatteringPhysical&);
  G4CMPIVScatteringPhysical& operator=(const G4CMPIVScatteringPhysical& right);
};

#endif

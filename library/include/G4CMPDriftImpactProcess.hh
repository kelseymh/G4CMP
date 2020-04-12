/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// 20200331  G4CMP-196: Added impact ionization process

#ifndef G4CMPDriftImpactProcess_h
#define G4CMPDriftImpactProcess_h 1

#include "G4CMPVDriftProcess.hh"


class G4CMPDriftImpactProcess : public G4CMPVDriftProcess {

public:
  G4CMPDriftImpactProcess(const G4String& name = "G4CMPChargeImpact");
  virtual ~G4CMPDriftImpactProcess();

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

protected:
  virtual G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*);

private:
  // No copying/moving
  G4CMPDriftImpactProcess(G4CMPDriftImpactProcess&);
  G4CMPDriftImpactProcess(G4CMPDriftImpactProcess&&);
  G4CMPDriftImpactProcess& operator=(const G4CMPDriftImpactProcess&);
  G4CMPDriftImpactProcess& operator=(const G4CMPDriftImpactProcess&&);
};

#endif	/* G4CMPDriftImpactProcess_h */

/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef G4CMPDriftRecombinationProcess_h
#define G4CMPDriftRecombinationProcess_h 1

#include "G4CMPVDriftProcess.hh"
#include "globals.hh"

class G4CMPDriftRecombinationProcess : public G4CMPVDriftProcess {
public:
  G4CMPDriftRecombinationProcess(const G4String& name = "G4CMPChargeRecombine",
                                 G4CMPProcessSubType type = fChargeRecombine);
  // No copying/moving
  G4CMPDriftRecombinationProcess(G4CMPDriftRecombinationProcess&) = delete;
  G4CMPDriftRecombinationProcess(G4CMPDriftRecombinationProcess&&) = delete;
  G4CMPDriftRecombinationProcess& operator=(const G4CMPDriftRecombinationProcess&) = delete;
  G4CMPDriftRecombinationProcess& operator=(const G4CMPDriftRecombinationProcess&&) = delete;

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&) override;

protected:
  virtual G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*)
                                   override;
};

#endif

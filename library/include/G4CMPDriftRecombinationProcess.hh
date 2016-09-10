/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef G4CMPDriftRecombinationProcess_h
#define G4CMPDriftRecombinationProcess_h 1

#include "G4CMPVDriftProcess.hh"
#include "G4VRestProcess.hh"
#include "globals.hh"

class G4CMPDriftRecombinationProcess : public G4VRestProcess, G4CMPProcessUtils {
public:
  G4CMPDriftRecombinationProcess(const G4String& name = "G4CMPChargeRecombine",
                                 G4CMPProcessSubType type = fChargeRecombine);

  virtual G4bool IsApplicable(const G4ParticleDefinition& aPD);
  virtual void StartTracking(G4Track* track);
  virtual void EndTracking();

  virtual G4double AtRestGetPhysicalInteractionLength(const G4Track& track,
                                                   G4ForceCondition* condition);

  virtual G4VParticleChange* AtRestDoIt(const G4Track&, const G4Step&);

protected:
  virtual G4double GetMeanLifeTime(const G4Track&, G4ForceCondition*);

private:
  // No copying/moving
  G4CMPDriftRecombinationProcess(G4CMPDriftRecombinationProcess&);
  G4CMPDriftRecombinationProcess(G4CMPDriftRecombinationProcess&&);
  G4CMPDriftRecombinationProcess& operator=(const G4CMPDriftRecombinationProcess&);
  G4CMPDriftRecombinationProcess& operator=(const G4CMPDriftRecombinationProcess&&);
};

#endif

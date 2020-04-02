/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// 20200331  G4CMP-195/196: Added impact ionization and trapping

#include "G4CMPDriftTrappingProcess.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4CMPEnergyPartition.hh"
#include "G4CMPSecondaryUtils.hh"
#include "G4CMPUtils.hh"
#include "G4LatticePhysical.hh"
#include "G4RandomDirection.hh"
#include "G4Track.hh"
#include <vector>


// Constructor and destructor

G4CMPDriftTrappingProcess::
G4CMPDriftTrappingProcess(const G4String &name)
  : G4CMPVDriftProcess(name, fTrapping) {}

G4CMPDriftTrappingProcess::~G4CMPDriftTrappingProcess() {
  delete partitioner;
}


// Process actions

G4double 
G4CMPDriftTrappingProcess::GetMeanFreePath(const G4Track&, G4double,
						G4ForceCondition* cond) {
  // CS TODO: How to access impactLength from ConfigManager?
  *cond = Forced;
  return DBL_MAX;
}

G4VParticleChange* 
G4CMPDriftTrappingProcess::PostStepDoIt(const G4Track& aTrack,
					     const G4Step& aStep) {
  aParticleChange.Initialize(aTrack);

  if (verboseLevel > 1) {
    G4cout << "G4CMPDriftTrappingProcess::PostStepDoIt: "
           << aTrack.GetDefinition()->GetParticleName()
           << " trapped by an impurity."
           << G4endl;
  }

  aParticleChange.ProposeTrackStatus(fStopAndKill);

  ClearNumberOfInteractionLengthLeft();		// All processes should do this!
  return &aParticleChange;
}

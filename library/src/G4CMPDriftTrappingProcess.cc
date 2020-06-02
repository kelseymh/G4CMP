/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// 20200331  C. Stanford G4CMP-195: Added charge trapping
// 20200411  M. Kelsey: Clean up undefined functions, retrieve preset MFP
// 20200331  G4CMP-195:  Add chage trapping MFP
// 20200504  G4CMP-195:  Reduce length of charge-trapping parameter names;
//		provide static function for MFP access; remove unnecessary
//		#includes.

#include "G4CMPDriftTrappingProcess.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPEnergyPartition.hh"
#include "G4CMPUtils.hh"
#include "G4Track.hh"


// Constructor and destructor

G4CMPDriftTrappingProcess::
G4CMPDriftTrappingProcess(const G4String &name)
  : G4CMPVDriftProcess(name, fChargeTrapping), partitioner(0) {
  // If needed, create G4CMPEnergyPartition here
}

G4CMPDriftTrappingProcess::~G4CMPDriftTrappingProcess() {
  delete partitioner;
}


// Evaluate MFP for current particle type

G4double 
G4CMPDriftTrappingProcess::GetMeanFreePath(const G4ParticleDefinition* pd) {
  return (G4CMP::IsElectron(pd) ? G4CMPConfigManager::GetETrappingMFP()
	  : G4CMP::IsHole(pd)   ? G4CMPConfigManager::GetHTrappingMFP()
	  : DBL_MAX);
}

G4double G4CMPDriftTrappingProcess::GetMeanFreePath(const G4Track&, G4double,
						    G4ForceCondition*) {
  return GetMeanFreePath(GetCurrentParticle());
}

// Process actions

G4VParticleChange* 
G4CMPDriftTrappingProcess::PostStepDoIt(const G4Track& aTrack,
					const G4Step& /*aStep*/) {
  aParticleChange.Initialize(aTrack);

  if (verboseLevel > 1) {
    G4cout << "G4CMPDriftTrappingProcess::PostStepDoIt: "
           << aTrack.GetDefinition()->GetParticleName()
           << " trapped by an impurity.  No energy released."
           << G4endl;
  }

  // NOTE: If trap depth allows for energy release, use partitioner here
  aParticleChange.ProposeTrackStatus(fStopAndKill);

  ClearNumberOfInteractionLengthLeft();		// All processes should do this!
  return &aParticleChange;
}

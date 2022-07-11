/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// 20170620  M. Kelsey -- Follow interface changes in G4CMPSecondaryUtils
// 20170802  M. Kelsey -- Replace phonon production with G4CMPEnergyPartition
// 20180827  M. Kelsey -- Prevent partitioner from recomputing sampling factors
// 20210328  Modify above; compute direct-phonon sampling factor here

#include "G4CMPDriftRecombinationProcess.hh"
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

G4CMPDriftRecombinationProcess::
G4CMPDriftRecombinationProcess(const G4String &name, G4CMPProcessSubType type)
  : G4CMPVDriftProcess(name, type), partitioner(new G4CMPEnergyPartition) {
  partitioner->UseDownsampling(false);		// Apply preset scaling factors
}

G4CMPDriftRecombinationProcess::~G4CMPDriftRecombinationProcess() {
  delete partitioner;
}


// Process actions

G4double 
G4CMPDriftRecombinationProcess::GetMeanFreePath(const G4Track&, G4double,
						G4ForceCondition* cond) {
  *cond = Forced;
  return DBL_MAX;
}

G4VParticleChange* 
G4CMPDriftRecombinationProcess::PostStepDoIt(const G4Track& aTrack,
					     const G4Step& aStep) {
  aParticleChange.Initialize(aTrack);

  // If the particle has not come to rest, do nothing
  if (aTrack.GetTrackStatus() != fStopButAlive) {
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }

  if (verboseLevel > 1) {
    G4cout << "G4CMPDriftRecombinationProcess::PostStepDoIt: "
           << aTrack.GetDefinition()->GetParticleName()
           << " reabsorbed by the lattice."
           << G4endl;
  }

  *(G4CMPProcessUtils*)partitioner = *(G4CMPProcessUtils*)this;
  partitioner->SetVerboseLevel(verboseLevel);
  partitioner->UseVolume(aTrack.GetVolume());

  // FIXME: Each charge carrier is independent, so it only gives back 0.5 times
  // the band gap. Really electrons and holes should recombine, killing both
  // tracks and giving back the band gap. Maybe there is a better way?
  G4double ePot = 0.5 * theLattice->GetBandGapEnergy();

  partitioner->ComputePhononSampling(ePot);
  partitioner->DoPartition(0., ePot);
  partitioner->GetSecondaries(&aParticleChange);

  if (aParticleChange.GetNumberOfSecondaries() == 0) {	// Record energy release
    aParticleChange.ProposeNonIonizingEnergyDeposit(ePot);
  }

  aParticleChange.ProposeTrackStatus(fStopAndKill);

  ClearNumberOfInteractionLengthLeft();		// All processes should do this!
  return &aParticleChange;
}

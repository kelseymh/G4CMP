/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPSecondaryProduction.hh
/// \brief Definition of the G4CMPSecondaryProduction process class.  This
///	class will be used to extend the existing Geant4 ionization
///	(and possibly other) processes to generate phonons and charge
///	carriers as secondaries.
//
// $Id$
//
// 20150306  Michael Kelsey
// 20160825  Replace implementation with use of G4CMPEnergyPartition
// 20191007  All normal G4 tracks should be used, not just charged.
// 20200222  Enable collection of EnergyPartition summary data.
// 20201207  Suspend parent track so that secondaries get processed first.
// 20210203  G4CMP-241 : Process must run after PostStepDoIt, not AlongStep.
// 20210303  G4CMP-243 : Consolidate nearby steps into one effective hit.
// 20210318  G4CMP-245 : Enforce clearance from crystal surfaces.
// 20210513  G4CMP-258 : Ensure that track weights are used with secondaries.
// 20210608  G4CMP-260 : Improve logic to collect steps and process hits in
//	       cases where a step doesn't have energy deposited.
// 20210610  G4CMP-262 : Handle step accumulation including track suspension,
//	       by keeping a map of accumulators by track ID
// 20220216  G4CMP-290 : Only spread secondaries along trajectory for charged
//	       tracks; neutrals get everything at endpoint.
// 20220815  G4CMP-308 : Factor step-accumulation procedures to HitMerging.
// 20220828  Call HitMerging::ProcessEvent() to ensure event ID is set.

#include "G4CMPSecondaryProduction.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPHitMerging.hh"
#include "G4CMPProcessSubType.hh"
#include "G4CMPUtils.hh"
#include "G4LatticeManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessType.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VParticleChange.hh"


// Constructor and destructor

G4CMPSecondaryProduction::G4CMPSecondaryProduction()
  : G4CMPVProcess("G4CMPSecondaryProduction", fSecondaryProduction),
    mergeHits(new G4CMPHitMerging), secondariesFirst(true) {;}

G4CMPSecondaryProduction::~G4CMPSecondaryProduction() {
  delete mergeHits;
}


// Applies to all charged, non-resonance particles except the drift charges

G4bool G4CMPSecondaryProduction::IsApplicable(const G4ParticleDefinition& pd) {
  return (!pd.IsShortLived() && !G4CMP::IsPhonon(pd) &&
	  !G4CMP::IsChargeCarrier(pd) );
}


// Override G4CMPProcessUtils for normal tracks outside lattice volumes

void G4CMPSecondaryProduction::LoadDataForTrack(const G4Track* track) {
  if (verboseLevel>1)
    G4cout << "G4CMPSecondaryProduction::LoadDataForTrack" << G4endl;

  SetCurrentTrack(track);

  // Skip further configuration if not active volume
  if (!G4LatticeManager::GetLatticeManager()->HasLattice(GetCurrentVolume())) {
    theLattice = 0;
    return;
  }

  SetLattice(track);

  // Configure hit merging
  mergeHits->SetVerboseLevel(verboseLevel);
  mergeHits->ProcessEvent();
  mergeHits->LoadDataForTrack(track);
}


// Use previously computed energy loss to generate secondaries

G4VParticleChange* 
G4CMPSecondaryProduction::PostStepDoIt(const G4Track& track,
				       const G4Step& step) {
  aParticleChange.Initialize(track); 

  // Only apply to tracks while they are in lattice-configured volumes
  LoadDataForTrack(&track);
  if (!theLattice) return &aParticleChange;

  if (verboseLevel) {
    G4cout << GetProcessName() << "::PostStepDoIt track " << &track
	   << " step " << &step << G4endl;
  }

  if (mergeHits->ProcessStep(step))
    mergeHits->FillOutput(&aParticleChange);

  // If requested (default), process new secondaries immediately
  if (aParticleChange.GetNumberOfSecondaries() > 0 &&
      secondariesFirst && track.GetTrackStatus() == fAlive)
    aParticleChange.ProposeTrackStatus(fSuspend);
  
  // NOTE:  This process does NOT change the track's momentum or energy
  return &aParticleChange;
}


// Process must be applied to all tracks at the end of their step

G4double 
G4CMPSecondaryProduction::GetMeanFreePath(const G4Track&, G4double,
					  G4ForceCondition* condition) {
  *condition = StronglyForced;
  return DBL_MAX;
}

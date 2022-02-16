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

#include "G4CMPSecondaryProduction.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4CMPEnergyPartition.hh"
#include "G4CMPGeometryUtils.hh"
#include "G4CMPProcessSubType.hh"
#include "G4CMPStepAccumulator.hh"
#include "G4CMPUtils.hh"
#include "G4Event.hh"
#include "G4IonisParamMat.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessType.hh"
#include "G4RandomDirection.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VParticleChange.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"
#include "Randomize.hh"
#include <algorithm>
#include <vector>


// Constructor and destructor

G4CMPSecondaryProduction::G4CMPSecondaryProduction()
  : G4CMPVProcess("G4CMPSecondaryProduction", fSecondaryProduction),
    partitioner(new G4CMPEnergyPartition),
    secondariesFirst(true), combiningStepLength(0.), accumulator(0),
    currentEventID(-1) {
  partitioner->FillSummary(true);	// Collect partition summary data
}

G4CMPSecondaryProduction::~G4CMPSecondaryProduction() {
  trackAccum.clear();
  delete partitioner;
}


// Applies to all charged, non-resonance particles except the drift charges

G4bool G4CMPSecondaryProduction::IsApplicable(const G4ParticleDefinition& pd) {
  return (!pd.IsShortLived() && !G4CMP::IsPhonon(pd) &&
	  !G4CMP::IsChargeCarrier(pd) );
}


// Overload G4CMPProcessUtils function to fill energy parameters

void G4CMPSecondaryProduction::LoadDataForTrack(const G4Track* track) {
  if (verboseLevel>1)
    G4cout << "G4CMPSecondaryProduction::LoadDataForTrack" << G4endl;

  // If new event, clear out any existing accumulators
  G4int thisEvent = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  if (thisEvent != currentEventID) {
    if (verboseLevel>1) G4cout << " New event: clearing accumulators" << G4endl;
    trackAccum.clear();
    currentEventID = thisEvent;
  }

  SetCurrentTrack(track);

  // Skip further configuration if not active volume
  if (!G4LatticeManager::GetLatticeManager()->HasLattice(GetCurrentVolume())) {
    theLattice = 0;
    return;
  }

  SetLattice(track);

  // Direct step accumulator to work with current track
  accumulator = &trackAccum[track->GetTrackID()];

  // Set up energy partitioning to work with current track and volume
  *(G4CMPProcessUtils*)partitioner = *(G4CMPProcessUtils*)this;
  partitioner->UseVolume(GetCurrentVolume());
  partitioner->SetVerboseLevel(verboseLevel);
}


// Use previously computed energy loss to generate secondaries

G4VParticleChange* 
G4CMPSecondaryProduction::PostStepDoIt(const G4Track& track,
				       const G4Step& stepData) {
  aParticleChange.Initialize(track); 
  LoadDataForTrack(&track);

  // Only apply to tracks while they are in lattice-configured volumes
  if (!theLattice) return &aParticleChange;

  if (verboseLevel) G4cout << GetProcessName() << "::PostStepDoIt" << G4endl;

  // Get configuration for how to merge steps
  combiningStepLength = G4CMPConfigManager::GetComboStepLength();

  // Check if current hit should be accumulated
  G4bool usedStep = DoAddStep(stepData);
  if (usedStep) {
    if (verboseLevel>1) {
      G4cout << " accumulating track " << track.GetTrackID()
	     << " step " << track.GetCurrentStepNumber()
	     << " @ " << stepData.GetPostStepPoint()->GetPosition()
	     << " Edep " << stepData.GetTotalEnergyDeposit()/eV << " eV"
	     << " Eniel " << stepData.GetNonIonizingEnergyDeposit()/eV << " eV"
	     << G4endl;
    }

    // Accumulate steps with non-zero energy deposit
    if (HasEnergy(stepData)) accumulator->Add(stepData);
  }

  // Check if effective hit should be generated
  G4bool generated = DoSecondaries(stepData);
  if (generated) {
    AddSecondaries();
    accumulator->Clear();
  }

  // If step wasn't added above, add it here for next time
  if (!usedStep && HasEnergy(stepData)) accumulator->Add(stepData);

  // If requested (default), process new secondaries immediately
  if (generated && secondariesFirst && track.GetTrackStatus() == fAlive)
    aParticleChange.ProposeTrackStatus(fSuspend);

  // NOTE:  This process does NOT change the track's momentum or energy
  return &aParticleChange;
}


// Decide if current step should be added to the accumulator

G4bool G4CMPSecondaryProduction::DoAddStep(const G4Step& stepData) const {
  G4StepStatus  sStatus = stepData.GetPostStepPoint()->GetStepStatus();
  G4TrackStatus tStatus = stepData.GetTrack()->GetTrackStatus();

  if (verboseLevel>1) {
    G4cout << " DoAddStep:"
	   << " nsteps==0 ? " << (accumulator->nsteps==0)
	   << "\n stepLen ? " << (stepData.GetStepLength()<combiningStepLength)
	   << "\n boundary ? " << (sStatus == fGeomBoundary ||
				   sStatus == fWorldBoundary)
	   << "\n stopped ? " << (tStatus != fAlive)
	   << " : tStatus " << tStatus << G4endl;
  }

  return (accumulator->nsteps == 0 ||	// First step always goes in
	  stepData.GetStepLength() < combiningStepLength ||
	  sStatus == fGeomBoundary || sStatus == fWorldBoundary ||
	  tStatus != fAlive);		// Stopping tracks must be caught
}

// Check if step includes any energy deposit

G4bool G4CMPSecondaryProduction::HasEnergy(const G4Step& stepData) const {
  G4double Edep  = stepData.GetTotalEnergyDeposit();
  G4double Eniel = stepData.GetNonIonizingEnergyDeposit();

  if (verboseLevel>1) {
    G4cout << " HasEnergy: Edep " << Edep/eV << " eV"
	   << " Eniel " << Eniel/eV << " eV" << G4endl;
  }

  return (Edep>0. || Eniel>0.);
}

// Decide if accumulator should be converted to secondaries

G4bool G4CMPSecondaryProduction::DoSecondaries(const G4Step& stepData) const {
  G4StepStatus  sStatus = stepData.GetPostStepPoint()->GetStepStatus();
  G4TrackStatus tStatus = stepData.GetTrack()->GetTrackStatus();

  if (verboseLevel>1) {
    G4cout << " DoSecondaries:"
	   << " nsteps>0 ? " << (accumulator->nsteps>0)
	   << "\n stepLen ? " << (stepData.GetStepLength()>=combiningStepLength)
	   << "\n boundary ? " << (sStatus == fGeomBoundary ||
				   sStatus == fWorldBoundary)
	   << "\n stopped ? " << (tStatus != fAlive)
	   << " : tStatus " << tStatus << G4endl;
  }

  return (accumulator->nsteps > 0 &&	// Don't process empty accumulator
	  (stepData.GetStepLength() >= combiningStepLength ||
	   sStatus == fGeomBoundary || sStatus == fWorldBoundary ||
	   tStatus != fAlive)		// Stopping tracks must be caught
	  );
}


// Use energy loss to generate phonons and charge carriers along path

void G4CMPSecondaryProduction::AddSecondaries() {
  G4double eTotal = accumulator->Edep;
  G4double eNIEL  = accumulator->Eniel;

  if (eTotal <= 0. && eNIEL <= 0.) return;	// Avoid unncessary work

  if (verboseLevel) {
    G4cout << " AddSecondaries from " << accumulator->nsteps << " steps "
	   << eTotal/eV << " eV" << " (" << eNIEL/eV << " NIEL)" << G4endl;
  }

  // Configure energy partitioning for EM, nuclear, or pre-determined energy
  G4int ptype = accumulator->pd->GetPDGEncoding();
  partitioner->DoPartition(ptype, eTotal, eNIEL);
  partitioner->GetSecondaries(theSecs);

  size_t nsec = theSecs.size();
  if (nsec == 0) {				// Avoid unnecessary work
    if (verboseLevel>1) G4cout << " No secondaries generated." << G4endl;
    return;
  }

  // Charged particles have energy spread along trajectory, from dE/dx
  if (accumulator->pd->GetPDGCharge() != 0) {
    if (verboseLevel>2)
      G4cout << " Charged track; spreading dE/dx along trajectory" << G4endl;

    std::random_shuffle(theSecs.begin(), theSecs.end(), RandomIndex);
    GeneratePositions(nsec, accumulator->start, accumulator->end);
  } else {
    GeneratePositions(nsec, accumulator->end, accumulator->end);
  }

  if (verboseLevel>1) G4cout << " Adding " << nsec << " secondaries" << G4endl;
  aParticleChange.SetNumberOfSecondaries(nsec);
  aParticleChange.SetSecondaryWeightByProcess(true);

  // Distribute generated particles along (straight line) trajectory
  for (size_t i=0; i<theSecs.size(); i++) {
    theSecs[i]->SetPosition(SurfaceClearance(posSecs[i]));
    aParticleChange.AddSecondary(theSecs[i]);

    if (verboseLevel>2) {
      G4Track* aSec = theSecs[i];
      G4cout << " secondary " << i << " : "
	     << aSec->GetParticleDefinition()->GetParticleName()
	     << " " << aSec->GetKineticEnergy()/eV << " eV "
	     << " along " << aSec->GetMomentumDirection()
	     << " @ " << aSec->GetPosition()
	     << " (wt " << aSec->GetWeight() << ")"
	     << G4endl;
    }
  }
}


// Generate intermediate points along step trajectory (straight line!)
// NOTE:  For MSC type deposition, these points ought to be a random walk

void G4CMPSecondaryProduction::GeneratePositions(size_t nsec,
						 const G4ThreeVector& start,
						 const G4ThreeVector& end) {
  if (verboseLevel>1) G4cout << " GeneratePositions " << nsec << G4endl;

  posSecs.clear();

  // If everything happens at a point, just fill the position vector
  if (start == end) {
    posSecs.resize(nsec, end);
    return;
  }

  // Get average distance between secondaries along (straight) trajectory
  G4ThreeVector traj = end - start;
  G4ThreeVector tdir = traj.unit();

  G4double length = traj.mag();
  G4double dl = length / G4double(nsec);
  G4double sigl = dl/6.;

  if (verboseLevel>1) {
    G4cout << " Choosing positions along " << length/mm << " mm " << tdir
	   << ": steps " << dl << " +- " << sigl << " mm" << G4endl;
  }

  posSecs.reserve(nsec);	// Avoid re-allocation memory churn

  G4double substep = 0.;
  G4ThreeVector lastPos = start;
  for (size_t i=0; i<nsec; i++) {
    substep = G4RandGauss::shoot(dl, sigl);
    lastPos += substep*tdir;
    posSecs.push_back(lastPos);
  }
}


// Ensure that secondaries all start inside volume

G4ThreeVector 
G4CMPSecondaryProduction::SurfaceClearance(const G4ThreeVector& pos) {
  return G4CMP::ApplySurfaceClearance(GetCurrentTouchable(), pos);
}


// Process must be applied to all tracks at the end of their step

G4double 
G4CMPSecondaryProduction::GetMeanFreePath(const G4Track&, G4double,
					  G4ForceCondition* condition) {
  *condition = StronglyForced;
  return DBL_MAX;
}


// Generate random index for shuffling secondaries

size_t G4CMPSecondaryProduction::RandomIndex(size_t n) {
  return (size_t)(n*G4UniformRand());
}

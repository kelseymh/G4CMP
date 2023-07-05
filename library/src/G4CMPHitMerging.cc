/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPHitMerging.hh
/// \brief Definition of the G4CMPHitMerging factory class.  This class
///	is used by G4CMPSecondaryProduction, and may be used directly by
///     user applications, to consolidate adjacent "hits" (energy deposits)
///     along a track, for more efficient downsampling.
//
// $Id$
//
// 20220815  Michael Kelsey -- Extracted from G4CMPSecondaryProduction
// 20220821  G4CMP-308 -- Use new G4CMPStepInfo container instead of G4Step
// 20220828  Pass event ID through to accumulator; improve debugging output

#include "G4CMPHitMerging.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPEnergyPartition.hh"
#include "G4CMPGeometryUtils.hh"
#include "G4CMPStepAccumulator.hh"
#include "G4CMPUtils.hh"
#include "G4Event.hh"
#include "G4Exception.hh"
#include "G4ExceptionSeverity.hh"
#include "G4LatticeManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VParticleChange.hh"
#include "G4VPhysicalVolume.hh"
#include "Randomize.hh"
#include <algorithm>
#include <vector>


// Constructor and destructor

G4CMPHitMerging::G4CMPHitMerging()
  : G4CMPProcessUtils(), verboseLevel(G4CMPConfigManager::GetVerboseLevel()),
    combiningStepLength(G4CMPConfigManager::GetComboStepLength()),
    accumulator(0), currentEventID(-1),
    partitioner(new G4CMPEnergyPartition) {
  partitioner->FillSummary(true);	// Collect partition summary data
}

G4CMPHitMerging::~G4CMPHitMerging() {
  trackAccum.clear();
  delete partitioner;
}


// Overload G4CMPProcessUtils function to initialize energy partitioner

void G4CMPHitMerging::LoadDataForTrack(const G4Track* track) {
  if (verboseLevel>1)
    G4cout << "G4CMPHitMerging::LoadDataForTrack" << G4endl;

  // If new event, clear out any existing accumulators
  ProcessEvent(G4RunManager::GetRunManager()->GetCurrentEvent());

  SetCurrentTrack(track);

  // Skip further configuration if not active volume
  if (!G4LatticeManager::GetLatticeManager()->HasLattice(GetCurrentVolume())) {
    theLattice = 0;
    return;
  }

  SetLattice(track);
}


// For primary generators, event must be passed in, not available from RM

void G4CMPHitMerging::ProcessEvent(const G4Event* currentEvent) {
  // If no event available, abort job with explanatory message
  if (!currentEvent)
    currentEvent = G4RunManager::GetRunManager()->GetCurrentEvent();

  if (!currentEvent) {
    currentEventID = -1;

    G4ExceptionDescription msg;
    msg << "No event available to start processing.  If called from\n"
	<< "primary generator action, ensure that ProcessEvent() is\n"
	<< "called with the local G4Event*.";
    G4Exception("G4CMPHitMerging::ProcessEvent", "Merging001",
		FatalException, msg);
    return;
  }
		
  // If new event, clear out any existing accumulators
  G4int thisEvent = currentEvent->GetEventID();
  if (thisEvent != currentEventID) {
    if (verboseLevel>1) G4cout << " New event: clearing accumulators" << G4endl;
    trackAccum.clear();
    currentEventID = thisEvent;
  }

  // Get configuration for how to merge steps
  combiningStepLength = G4CMPConfigManager::GetComboStepLength();
  if (verboseLevel>1) {
    G4cout << " combining steps within " << combiningStepLength << " mm"
	   << G4endl;
  }
}


// Use previously computed energy loss to generate secondaries

G4bool G4CMPHitMerging::ProcessStep(const G4Step& step) {
  if (verboseLevel>1)
    G4cout << "G4CMPHitMerging::ProcessStep(const G4Step&) " << &step << G4endl;

  return ProcessStep(G4CMPStepInfo(step));
}

G4bool G4CMPHitMerging::ProcessStep(const G4Step* step) {
  if (verboseLevel>1)
    G4cout << "G4CMPHitMerging::ProcessStep(const G4Step*) " << step << G4endl;

  return (step && ProcessStep(*step));
}

G4bool G4CMPHitMerging::ProcessStep(const G4CMPStepInfo& stepData) {
  // Only apply to tracks while they are in lattice-configured volumes
  readyForOutput = false;
  if (!theLattice) return false;

  if (verboseLevel) G4cout << "G4CMPHitMerging::ProcessStep" << G4endl;

  // Direct step accumulator to work with current track and event
  accumulator = &trackAccum[stepData.trackID];	// Creates new if needed
  accumulator->ProcessEvent(currentEventID);

  // Set up energy partitioning to work with current track and volume
  *(G4CMPProcessUtils*)partitioner = *(G4CMPProcessUtils*)this;
  partitioner->SetVerboseLevel(verboseLevel);
  partitioner->UseVolume(GetCurrentVolume());

  // Get configuration for how to merge steps
  combiningStepLength = G4CMPConfigManager::GetComboStepLength();

  // Check if current hit should be accumulated
  G4bool usedStep = DoAddStep(stepData);
  if (usedStep) {
    if (verboseLevel>1) {
      G4cout << " accumulating track " << stepData.trackID
	     << " step " << stepData.stepID
	     << " @ " << stepData.end
	     << " Edep " << stepData.Edep/eV << " eV"
	     << " Eniel " << stepData.Eniel/eV << " eV"
	     << G4endl;
    }

    // Accumulate steps with non-zero energy deposit
    if (HasEnergy(stepData)) accumulator->Add(stepData);
  }

  // Check if effective hit should be generated
  readyForOutput = ReadyForOutput(stepData);
  if (readyForOutput) {
    PrepareOutput();
    accumulator->Clear();
  }

  // If step wasn't added above, add it here for next time
  if (!usedStep && HasEnergy(stepData)) accumulator->Add(stepData);

  return readyForOutput;
}


// Decide if current step should be added to the accumulator

G4bool G4CMPHitMerging::DoAddStep(const G4CMPStepInfo& stepData) const {
  if (verboseLevel>2) G4cout << " DoAddStep " << stepData << G4endl;

  G4double distance = (stepData.end - accumulator->end).mag();

  if (verboseLevel>1) {
    G4cout << " DoAddStep returns OR of the following: "
	   << "\n nsteps=0 ? " << (accumulator->nsteps==0)
	   << "\n dist<max ? " << (distance < combiningStepLength)
	   << "\n boundary ? " << (stepData.sStatus == fGeomBoundary ||
				   stepData.sStatus == fWorldBoundary)
	   << "\n stopped  ? " << (stepData.tStatus != fAlive)
	   << " : tStatus " << stepData.tStatus << G4endl;
  }

  return (accumulator->nsteps == 0 ||	// First step always goes in
	  distance < combiningStepLength ||
	  stepData.sStatus == fGeomBoundary ||
	  stepData.sStatus == fWorldBoundary ||
	  stepData.tStatus != fAlive);		// Stopping tracks must be caught
}

// Check if step includes any energy deposit

G4bool G4CMPHitMerging::HasEnergy(const G4CMPStepInfo& stepData) const {
  if (verboseLevel>1) {
    G4cout << " HasEnergy: Edep " << stepData.Edep/eV << " eV"
	   << " Eniel " << stepData.Eniel/eV << " eV" << G4endl;
  }

  return (stepData.Edep>0. || stepData.Eniel>0.);
}

// Decide if accumulator should be converted to secondaries

G4bool G4CMPHitMerging::ReadyForOutput(const G4CMPStepInfo& stepData) const {
  if (verboseLevel>2) G4cout << " ReadyForOutput " << stepData << G4endl;

  G4double distance = (stepData.end - accumulator->end).mag();

  if (verboseLevel>1) {
    G4cout << " ReadyForOutput return nsteps and OR of everything else:"
	   << "\n nsteps>0 ? " << (accumulator->nsteps>0)
	   << "\n dist>max ? " << (distance >= combiningStepLength)
	   << "\n boundary ? " << (stepData.sStatus == fGeomBoundary ||
				   stepData.sStatus == fWorldBoundary)
	   << "\n stopped  ? " << (stepData.tStatus != fAlive)
	   << " : tStatus " << stepData.tStatus << G4endl;
  }

  return (accumulator->nsteps > 0 &&	// Don't process empty accumulator
	  (distance >= combiningStepLength ||
	   stepData.sStatus == fGeomBoundary ||
	   stepData.sStatus == fWorldBoundary ||
	   stepData.tStatus != fAlive)		// Stopping tracks must be caught
	  );
}


// Use energy loss to generate phonons and charge carriers along path

void G4CMPHitMerging::PrepareOutput() {
  if (!readyForOutput) return;			// Not ready to do this
  if (accumulator->Edep <= 0. && accumulator->Eniel <= 0.) return;

  if (verboseLevel) {
    G4cout << "G4CMPHitMerging::PrepareOutput" << G4endl
	   << *accumulator << G4endl;
  }

  // Process recorded energy deposit(s) into phonons and charge carriers
  partitioner->DoPartition(accumulator);

  size_t nsec = partitioner->GetNumberOfTracks();
  if (nsec == 0) {				// Avoid unnecessary work
    if (verboseLevel>1) G4cout << " No secondaries generated." << G4endl;
    return;
  }

  // Charged particles have energy spread along trajectory, from dE/dx
  if (accumulator->pd->GetPDGCharge() != 0) {
    if (verboseLevel>2)
      G4cout << " Charged track; spreading dE/dx along trajectory" << G4endl;

    GeneratePositions(nsec, accumulator->start, accumulator->end);
  } else {
    GeneratePositions(1, accumulator->end, accumulator->end);
  }
}

// Populate particleChange with secondary tracks

void G4CMPHitMerging::FillOutput(G4VParticleChange* aParticleChange) {
  if (!readyForOutput) return;			// Nothing to be done

  partitioner->GetSecondaries(theSecs);
  size_t nsec = theSecs.size();
  if (nsec == 0) return;			// Nothing to be done

  if (verboseLevel>1) G4cout << " Adding " << nsec << " secondaries" << G4endl;

  aParticleChange->SetNumberOfSecondaries(nsec);
  aParticleChange->SetSecondaryWeightByProcess(true);

  // Distribute generated particles along positions within trajectory
  size_t npos = posSecs.size();
  for (size_t i=0; i<nsec; i++) {
    theSecs[i]->SetPosition(posSecs[std::min(npos-1,i)]);
    aParticleChange->AddSecondary(theSecs[i]);

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
  }	// for (i<nsec
}

// Populate event with primary tracks and vertices

void G4CMPHitMerging::FillOutput(G4Event* primaryEvent, G4double time) {
  if (!readyForOutput) return;			// Nothing to be done

  size_t nprim = partitioner->GetNumberOfTracks();
  if (nprim == 0) return;			// Nothing to be done

  if (verboseLevel>1) G4cout << " Adding " << nprim << " primaries" << G4endl;

  size_t npos = posSecs.size();

  // For a single point, let partitioner take care of event
  if (npos == 1) partitioner->GetPrimaries(primaryEvent, posSecs[0], time);
  else partitioner->GetPrimaries(primaryEvent, posSecs, time);
}


// Check for any non-empty accumulators, and generate primaries from them

void G4CMPHitMerging::FinishOutput(G4Event* primaryEvent) {
  if (trackAccum.empty()) return;		// Nothing to be done
  if (!primaryEvent) return;

  if (primaryEvent->GetEventID() != currentEventID) {
    G4ExceptionDescription msg;
    msg << "primaryEvent " << primaryEvent->GetEventID() << " does not"
	<< " match currentEventID " << currentEventID << ".  Accumulated"
	<< " hits may be lost.";
    G4Exception("G4CMPHitMerging::FinishOutput", "Merging002",
		JustWarning, msg);
    return;
  }

  if (verboseLevel) {
    G4cout << "G4CMPHitMerging::FinishOutput " << primaryEvent->GetEventID()
	   << G4endl;
  }

  // Loop over all registered accumulators and flush them to output
  for (const auto& trkAcc: trackAccum) {
    FlushAccumulator(trkAcc.first, primaryEvent);
  }
}

// Process specified accumulator into new primaries for event

void G4CMPHitMerging::FlushAccumulator(G4int trkID, G4Event* primaryEvent) {
  accumulator = &trackAccum[trkID];
  if (accumulator->nsteps == 0) return;		// Nothing to be done
  
  if (verboseLevel>1) {
    G4cout << "G4CMPHitMerging::FlushAccumulator track " << trkID << " with"
	   << accumulator->nsteps << " steps" << G4endl;
  }
    
  readyForOutput = true;
  PrepareOutput();
  FillOutput(primaryEvent, accumulator->time);

  accumulator->Clear();
}

// Generate intermediate points along step trajectory (straight line!)
// NOTE:  For MSC type deposition, these points ought to be a random walk

void G4CMPHitMerging::GeneratePositions(size_t nsec,
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
    lastPos += SurfaceClearance(substep*tdir);
    
    posSecs.push_back(lastPos);
  }
}


// Ensure that secondaries all start inside volume

G4ThreeVector 
G4CMPHitMerging::SurfaceClearance(const G4ThreeVector& pos) {
  return G4CMP::ApplySurfaceClearance(GetCurrentTouchable(), pos);
}

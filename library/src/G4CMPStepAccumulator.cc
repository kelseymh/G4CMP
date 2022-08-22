/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPStepAccumulator.hh
/// \brief Implementation of the G4CMPStepAccumulator container.  This class  
///	stores energy deposit information from track steps, to create
///     phonon and charge carrier secondaries.
//
// 20210303  Michael Kelsey
// 20210608  Reset trackID in Clear(), check for matching eventID, and report
//	       rollover errors for track or event changes.
// 20220216  Add "stepID" to provide full identification.  Add printout.
// 20220228  Add sanity check that only energy-deposit hits are accumulated.
// 20220821  G4CMP-308 -- Define step-info container to avoid needing G4Step

#include "globals.hh"
#include "G4CMPStepAccumulator.hh"
#include "G4Event.hh"
#include "G4ParticleDefinition.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include <iostream>


// Populate step-info container

G4CMPStepInfo::G4CMPStepInfo(const G4Step* step)
  : trackID(step->GetTrack()->GetTrackID()),
    stepID(step->GetTrack()->GetCurrentStepNumber()),
    pd(step->GetTrack()->GetDefinition()),
    Edep(step->GetTotalEnergyDeposit()),
    Eniel(step->GetNonIonizingEnergyDeposit()),
    start(step->GetPreStepPoint()->GetPosition()),
    end(step->GetPostStepPoint()->GetPosition()) {;}

G4CMPStepInfo::G4CMPStepInfo(const G4Step& step)
  : trackID(step.GetTrack()->GetTrackID()),
    stepID(step.GetTrack()->GetCurrentStepNumber()),
    pd(step.GetTrack()->GetDefinition()),
    Edep(step.GetTotalEnergyDeposit()),
    Eniel(step.GetNonIonizingEnergyDeposit()),
    start(step.GetPreStepPoint()->GetPosition()),
    end(step.GetPostStepPoint()->GetPosition()) {;}

// Duplicate information from specified other step

G4CMPStepInfo&
G4CMPStepInfo::operator=(const G4CMPStepInfo& step) {
  trackID = step.trackID;
  stepID = step.stepID;
  pd = step.pd;
  Edep = step.Edep;
  Eniel = step.Eniel;
  start = step.start;
  end = step.end;

  return *this;
}


// Extract relevant information from step

void G4CMPStepAccumulator::Add(const G4CMPStepInfo& step) {
  // If event ID has changed discard previous data and report error
  G4int thisEvt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  if (eventID >=0 && eventID != thisEvt) {
    G4cerr << "ERROR G4CMPStepAccumulator rolled over between events "
	   << eventID << " and " << thisEvt << G4endl
	   << " Energy " << Edep/eV << " eV,"
	   << " NIEL " << Eniel/eV << " eV"
	   << " lost from previous event." << G4endl;
    Clear();
  }

  // If track type or track ID has changed, discard previous data
  if (trackID >= 0 && (trackID != step.trackID || pd != step.pd)) {
    G4cerr << "ERROR G4CMPStepAccumulator rolled over between tracks "
	   << trackID << " and " << step.trackID << G4endl
           << " Energy " << Edep/eV << " eV,"
	   << " NIEL " << Eniel/eV << " eV"
           << " lost from previous track." << G4endl;
    Clear();
  }

  // Should never receive zero-energy step
  if (step.Edep <= 0. && step.Eniel <= 0.) {
    G4cerr << "ERROR G4CMPStepAccumulator received zero-energy step."
	   << G4endl;
    return;
  }

  // Initialize new accumulation
  if (nsteps == 0) {
    eventID = thisEvt;
    (G4CMPStepInfo)(*this) = step;
  }

  // Update step information
  nsteps++;
  stepID = step.stepID;		// Step ID will always point to last step

  // Compute energy-weighted centroid of step endpoints
  ((end *= Edep) += step.end*step.Edep) /= (Edep+step.Edep);
  // NOTE: Using accumulators above as lvalues to avoid creating temporaries
  // Equivalent to: end = (end*Edep + stepEnd*stepEdep)/(Edep+stepEdep);
  
  // Accumulate total energy of steps
  Edep  += step.Edep;
  Eniel += step.Eniel;
}


// Dump content for diagnostics

void G4CMPStepAccumulator::Print(std::ostream& os) const {
  os << "G4CMPStepAccumulator event " << eventID
     << " " << (pd?pd->GetParticleName():"")
     << " track " << trackID << "/" << stepID
     << "\n " << nsteps << " steps from " << start << " to " << end
     << "\n deposited " << Edep/eV << " eV, "
     << Eniel/eV << " eV NIEL"
     << std::endl;
}

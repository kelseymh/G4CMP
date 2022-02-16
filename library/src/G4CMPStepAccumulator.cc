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


// Reset accumulator for new track

void G4CMPStepAccumulator::Clear() {
  eventID = trackID = stepID = -1;
  nsteps = 0;
  pd = 0;
  Edep = Eniel = 0.;
  start.set(0,0,0);
  end.set(0,0,0);
}


// Extract relevant information from step

void G4CMPStepAccumulator::Add(const G4Step& step) {
  // If event ID has changed discard previous data and report error
  G4int thisEvt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  if (eventID >=0 && eventID != thisEvt) {
    G4cerr << "ERROR G4CMPStepAccumulator rolled over between events "
	   << eventID << " and " << thisEvt << G4endl
	   << " Energy " << Edep/eV << " eV, NIEL " << Eniel/eV << " eV"
	   << " lost from previous event." << G4endl;
    Clear();
  }
  
  // If track type or track ID has changed, discard previous data
  const G4Track* track = step.GetTrack();
  if (trackID >= 0 && (trackID != track->GetTrackID() ||
		       pd != track->GetParticleDefinition())) {
    G4cerr << "ERROR G4CMPStepAccumulator rolled over between tracks "
	   << trackID << " and " << track->GetTrackID() << G4endl
           << " Energy " << Edep/eV << " eV, NIEL " << Eniel/eV << " eV"
           << " lost from previous track." << G4endl;
    Clear();
  }
  
  if (nsteps == 0) {
    eventID = thisEvt;
    trackID = track->GetTrackID();
    pd = track->GetParticleDefinition();
    start = step.GetPreStepPoint()->GetPosition();
  }

  // Accumulate energy from current step
  nsteps++;
  stepID = track->GetCurrentStepNumber();
  Edep  += step.GetTotalEnergyDeposit();
  Eniel += step.GetNonIonizingEnergyDeposit();

  // Move endpoint to current step
  end = step.GetPostStepPoint()->GetPosition();
}


// Dump content for diagnostics

void G4CMPStepAccumulator::Print(std::ostream& os) const {
  os << "G4CMPStepAccumulator event " << eventID
     << " " << (pd?pd->GetParticleName():"")
     << " track " << trackID << "/" << stepID
     << "\n " << nsteps << " steps from " << start << " to " << end
     << "\n deposited " << Edep/eV << " eV, " << Eniel/eV << " eV NIEL"
     << std::endl;
}

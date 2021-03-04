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

#include "globals.hh"
#include "G4CMPStepAccumulator.hh"
#include "G4ParticleDefinition.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"


// Reset accumulator for new track

void G4CMPStepAccumulator::Clear() {
  nsteps = 0;
  pd = 0;
  Edep = Eniel = 0.;
  start.set(0,0,0);
  end.set(0,0,0);
}


// Extract relevant information from step

void G4CMPStepAccumulator::Add(const G4Step& step) {
  const G4Track* track = step.GetTrack();

  // If track type or track ID has changed, discard previous data
  if (trackID != track->GetTrackID() ||
      pd != track->GetParticleDefinition()) Clear();

  if (nsteps == 0) {
    trackID = track->GetTrackID();
    pd = track->GetParticleDefinition();
    start = step.GetPreStepPoint()->GetPosition();
  }

  // Accumulate energy from current step
  nsteps++;
  Edep  += step.GetTotalEnergyDeposit();
  Eniel += step.GetNonIonizingEnergyDeposit();

  // Move endpoint to current step
  end = step.GetPostStepPoint()->GetPosition();
}

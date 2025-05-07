/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPTrackLimiter.cc
/// \brief Implementation of the G4CMPTrackLimiter process, to kill tracks
///        falling below a minimum energy (set in G4CMPConfigManager).
//
// $Id$
//
// 20170822  M. Kelsey -- Add checking on current vs. original volume
// 20240506  G4CMP-371 -- Add flag to keep or discard below-minimum track energy
// 20250501  G4CMP-358 -- Identify and stop charge tracks stuck in field,
//	       using maxSteps configuration parameter.
// 20250506  Add local caches to compute cumulative flight distance, RMS

#include "G4CMPTrackLimiter.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPUtils.hh"
#include "G4ForceCondition.hh"
#include "G4ParticleChange.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include <limits.h>


// Only applies to G4CMP particles

G4bool G4CMPTrackLimiter::IsApplicable(const G4ParticleDefinition& pd) {
  return (G4CMP::IsPhonon(pd) || G4CMP::IsChargeCarrier(pd));
}


// Initialize flight-distance accumulators for new track

void G4CMPTrackLimiter::LoadDataForTrack(const G4Track* track) {
  G4CMPProcessUtils::LoadDataForTrack(track);

  flightAvg = flightAvg2 = lastFlight10k = lastRMS10k = 0.;
}


// Force killing if below cut

G4double G4CMPTrackLimiter::GetMeanFreePath(const G4Track&, G4double,
					    G4ForceCondition* condition) {
  *condition = StronglyForced;	// Ensures execution even with other Forced
  return DBL_MAX;
}

G4double G4CMPTrackLimiter::
PostStepGetPhysicalInteractionLength(const G4Track& trk, G4double sl,
				     G4ForceCondition* condition) {
  return GetMeanFreePath(trk, sl, condition);	// No GPIL handling needed
}

G4VParticleChange* G4CMPTrackLimiter::PostStepDoIt(const G4Track& track,
                                                    const G4Step& step) {
  aParticleChange.Initialize(track);

  if (verboseLevel>1) G4cout << GetProcessName() << "::PostStepDoIt" << G4endl;

  // Apply minimum energy cut to kill tracks with optional NIEL deposit
  if (BelowEnergyCut(track)) {
    if (verboseLevel>2) G4cout << " track below minimum energy." << G4endl;

    aParticleChange.ProposeTrackStatus(fStopAndKill);

    if (G4CMPConfigManager::RecordMinETracks())
      aParticleChange.ProposeNonIonizingEnergyDeposit(track.GetKineticEnergy());
  }

  // Ensure that track is still in original, valid volume
  if (EscapedFromVolume(step)) {
    G4Exception("G4CMPTrackLimiter", "Limit001", JustWarning,
		"Killing track escaped from original volume.");

    aParticleChange.SetNumberOfSecondaries(0);	// Don't launch bad tracks!
    aParticleChange.ProposeTrackStatus(fStopAndKill);
  }

  // Ensure track has not gotten stuck somewhere in mesh field
  if (ChargeStuck(track)) {
    std::stringstream msg;
    msg << "Stopping charged track stuck in mesh electric field @ "
	<< GetLocalPosition(track) << " local" << G4endl;
    G4Exception("G4CMPTrackLimiter", "Limit003", JustWarning,
		msg.str().c_str());

    aParticleChange.ProposeTrackStatus(fStopButAlive);
  }

  return &aParticleChange;
}


// Evaluate current track

G4bool G4CMPTrackLimiter::BelowEnergyCut(const G4Track& track) const {
  G4double ecut =
    (G4CMP::IsChargeCarrier(track) ? G4CMPConfigManager::GetMinChargeEnergy()
     : G4CMP::IsPhonon(track) ? G4CMPConfigManager::GetMinPhononEnergy() : -1.);

  return (track.GetKineticEnergy() < ecut);
}

G4bool G4CMPTrackLimiter::EscapedFromVolume(const G4Step& step) const {
  G4VPhysicalVolume* prePV  = step.GetPreStepPoint()->GetPhysicalVolume();
  G4VPhysicalVolume* postPV = step.GetPostStepPoint()->GetPhysicalVolume();

  if (verboseLevel>2) {
    G4cout << " prePV " << prePV->GetName()
	   << " postPV " << (postPV?postPV->GetName():"OutOfWorld")
	   << " status " << step.GetPostStepPoint()->GetStepStatus()
	   << G4endl;
  }

  // Track is NOT at a boundary, is stepping outside volume, or already escaped
  return ( (step.GetPostStepPoint()->GetStepStatus() != fGeomBoundary) &&
	   (postPV != GetCurrentVolume() || prePV != GetCurrentVolume())
	   );
}

// Note: non-const here to use accumulator caches

G4bool G4CMPTrackLimiter::ChargeStuck(const G4Track& track) {
  if (!IsChargeCarrier()) return false;		// Ignore phonons (for now?)

  // How long and how far has the track been travelling?
  const G4double maxSteps = G4CMPConfigManager::GetMaxChargeSteps();
  G4int nstep = track.GetCurrentStepNumber();

  G4double pathLen = track.GetTrackLength();
  G4double flightDist = (track.GetPosition()-track.GetVertexPosition()).mag();

  // Accumulate flight distance and and sum-of-squares averages
  flightAvg  = flightAvg + (flightDist-flightAvg)/nstep;
  flightAvg2 = flightAvg2 + (flightDist*flightDist - flightAvg2)/nstep;

  // Compute change in distance every 10,000 steps
  G4double fltChange = flightDist-lastFlight10k;
  if (nstep%10000 == 0) lastFlight10k = flightDist;

  // Compute change in RMS every 10,000 steps
  const G4double minRMS = 0.;			// Zero means never bad
  G4double RMS = sqrt(flightAvg2 - flightAvg*flightAvg);
  G4double RMSchange = RMS-lastRMS10k;
  if (nstep%10000 == 0) lastRMS10k = RMS;

  // Scattering makes the path length longer, but only a factor of a few
  const G4double maxScale = 20.;		// Not used
  G4double pathScale = pathLen / flightDist;

  if (verboseLevel>1) {
    G4cout << " after " << nstep << " steps, path " << pathLen
	   << " ~ " << pathScale << " x flight " << flightDist
	   << " (" << (pathScale>maxScale?">":"<") << maxScale << ")" << G4endl
	   << " changed by " << fltChange << " RMS " << RMS
	   << " changed by " << RMSchange << G4endl;
  }

  // 
  return ((maxSteps>0 && nstep>maxSteps) ||
	  (nstep>1000 && nstep%10000==0 && fabs(RMSchange) < minRMS));
}

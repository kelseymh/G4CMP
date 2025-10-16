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
// 20240506  G4CMP-371:  Add flag to keep or discard below-minimum track energy.
// 20250327  G4CMP-468:  Stop surface displacement reflections from "escaping."
// 20250413  G4CMP-468:  Move diagnostic outputs inside verbosity.
// 20250421  Add comparison of track position with volume.
// 20250501  G4CMP-358:  Identify and stop charge tracks stuck in field,
//	       using maxSteps configuration parameter.
// 20250506  Add local caches to compute cumulative flight distance, RMS
// 20250801  G4CMP-326:  Kill thermal phonons if finite temperature set.
// 20251015  G4CMP-516:  Add excessPath to ChargeStuck() boolean return.

#include "G4CMPTrackLimiter.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPGeometryUtils.hh"
#include "G4CMPUtils.hh"
#include "G4ForceCondition.hh"
#include "G4LatticePhysical.hh"
#include "G4Navigator.hh"
#include "G4ParticleChange.hh"
#include "G4PhysicalConstants.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4TransportationManager.hh"
#include "G4VSolid.hh"
#include <limits.h>
#include <sstream>


// Only applies to G4CMP particles

G4bool G4CMPTrackLimiter::IsApplicable(const G4ParticleDefinition& pd) {
  return (G4CMP::IsPhonon(pd) || G4CMP::IsChargeCarrier(pd));
}


// Initialize flight-distance accumulators for new track

void G4CMPTrackLimiter::LoadDataForTrack(const G4Track* track) {
  G4CMPProcessUtils::LoadDataForTrack(track);

  flightAvg = flightAvg2 = lastFlight = lastRMS = 0.;
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

  // Skip reflection zero-length steps
  if (step.GetStepLength() == 0.) return &aParticleChange;

  // Apply minimum energy cut to kill tracks with optional NIEL deposit
  if (BelowEnergyCut(track)) {
    if (verboseLevel>2) G4cout << " track below minimum energy." << G4endl;

    aParticleChange.ProposeTrackStatus(fStopAndKill);

    if (G4CMPConfigManager::RecordMinETracks())
      aParticleChange.ProposeNonIonizingEnergyDeposit(track.GetKineticEnergy());
  }

  // Ensure that track is still in original, valid volume
  if (EscapedFromVolume(step)) {
    std::stringstream msg;
    msg << "Killing track escaped from volume "
	<< GetCurrentVolume()->GetName() + ":"
	<< GetCurrentVolume()->GetCopyNo();
    G4Exception("G4CMPTrackLimiter", "Limit001", JustWarning,
		msg.str().c_str());

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

  // Check whether track position and volume are consistent
  if (InvalidPosition(track)) {
    std::stringstream msg;
    msg << "Killing track inconsistent position " << track.GetPosition()
	<< "\n vs. detector volume " << GetCurrentVolume()->GetName() + ":"
	<< GetCurrentVolume()->GetCopyNo();
    G4Exception("G4CMPTrackLimiter", "Limit002", JustWarning,
		msg.str().c_str());
    
    aParticleChange.SetNumberOfSecondaries(0);	// Don't launch bad tracks!
    aParticleChange.ProposeTrackStatus(fStopAndKill);
  }

  // Kill phonons consistent with thermal populations
  if (PhononIsThermal(track))
    aParticleChange.ProposeTrackStatus(fStopAndKill);
   
  return &aParticleChange;
}


// Evaluate current track

G4bool G4CMPTrackLimiter::BelowEnergyCut(const G4Track& track) const {
  G4double ecut =
    (G4CMP::IsChargeCarrier(track) ? G4CMPConfigManager::GetMinChargeEnergy()
     : G4CMP::IsPhonon(track) ? G4CMPConfigManager::GetMinPhononEnergy() : -1.);

  return (track.GetKineticEnergy() < ecut);
}

G4bool G4CMPTrackLimiter::InvalidPosition(const G4Track& track) const {
  G4VPhysicalVolume* trkVol = track.GetVolume();
  if (!trkVol) return false;

  const G4VTouchable* trkVT = track.GetTouchable();
  G4ThreeVector trkPos = track.GetPosition();
  if (verboseLevel>1) {
    G4cout << GetProcessName() << "::InvalidPosition()" << G4endl
	   << " trkVol " << trkVol->GetName() << " @ " << trkPos << G4endl;
  }

  G4CMP::RotateToLocalPosition(trkVT, trkPos);
  G4VSolid* solid = GetCurrentVolume()->GetLogicalVolume()->GetSolid();
  EInside isIn = solid->Inside(trkPos);
  if (verboseLevel>1) {
    const char* inName = (isIn==kInside ? "inside" : isIn==kOutside
			  ? "outside" : "surface");

    G4cout << " local " << trkPos << " is " << inName << " volume" << G4endl;
  }

  return (isIn == kOutside);
}

G4bool G4CMPTrackLimiter::EscapedFromVolume(const G4Step& step) const {
    G4StepPoint* preS = step.GetPreStepPoint();
    G4StepPoint* postS = step.GetPostStepPoint();

  G4VPhysicalVolume* prePV  = step.GetPreStepPoint()->GetPhysicalVolume();
  G4VPhysicalVolume* postPV = step.GetPostStepPoint()->GetPhysicalVolume();

  if (verboseLevel>1) {
    G4cout << GetProcessName() << "::EscapedFromVolume()" << G4endl
	   << " prePV " << prePV->GetName()
	   << " postPV " << (postPV?postPV->GetName():"OutOfWorld")
	   << " status " << postS->GetStepStatus()
	   << G4endl;

    if (verboseLevel>2) {
      const G4ThreeVector& prePt = preS->GetPosition();
      const G4ThreeVector& postPt = postS->GetPosition();

      G4cout << std::setprecision(std::numeric_limits<double>::max_digits10)
	     << "preStep status " << preS->GetStepStatus() << G4endl
	     << "preStep Pos    " << prePt << G4endl
	     << "postStep Pos   " << postPt << G4endl
	     << "stepPos dir    " << (postPt - prePt).unit() << G4endl
	     << "step Mom dir   " << postS->GetMomentumDirection() << G4endl
	     << "currPV Name    " << GetCurrentVolume()->GetName() << G4endl;

      G4VSolid* solid = GetCurrentVolume()->GetLogicalVolume()->GetSolid();
      EInside isIn = solid->Inside(GetLocalPosition(postS->GetPosition()));
      const char* inName = (isIn==kInside ? "inside" : isIn==kOutside
			    ? "outside" : "surface");
      G4cout << "Value for surface " << inName << G4endl;
    }
  }

  // Track is NOT at a boundary, is stepping outside volume, or already escaped
  G4bool escape =
    ((postS->GetStepStatus() != fGeomBoundary) &&
     (postPV != GetCurrentVolume() || prePV != GetCurrentVolume()));

  if (verboseLevel>1) G4cout << " escape? " << escape << G4endl;
  
  return escape;
}

G4bool G4CMPTrackLimiter::PhononIsThermal(const G4Track& track) const {
  if (!IsPhonon()) return false;	// Only phonons thermalize

  G4double temp = theLattice->GetTemperature();
  if (temp <= 0.) return false;

  if (verboseLevel>1)
    G4cout << GetProcessName() << "::PhononIsThermal()" << G4endl;

  if (verboseLevel>2) {
    G4cout << " phonon " << track.GetKineticEnergy()/eV << " eV"
	   << " vs. kT " << temp*k_Boltzmann/eV << " @ " << temp/kelvin << " K"
	   << G4endl;
  }

  G4bool isThermal = G4CMP::IsThermalized(track);
  if (verboseLevel>1) G4cout << " thermal? " << isThermal << G4endl;

  return isThermal;
}

// Note: non-const here to use accumulator caches

G4bool G4CMPTrackLimiter::ChargeStuck(const G4Track& track) {
  if (!IsChargeCarrier()) return false;		// Ignore phonons (for now?)

  // How long and how far has the track been travelling?
  const G4double maxSteps = G4CMPConfigManager::GetMaxChargeSteps();
  G4int nstep = track.GetCurrentStepNumber();

  G4double pathLen = track.GetTrackLength();
  const G4ThreeVector& pos = track.GetPosition();
  G4double flightDist = (pos-track.GetVertexPosition()).mag();
  G4double posShift = (nstep>1) ? (pos-lastPos).mag() : 0.;


  if (nstep%stepWindow == 1) {		// Start new window
    flightAvg = flightAvg2 = 0.;
    lastPos = pos;
  }

  // Accumulate flight distance and and sum-of-squares averages
  flightAvg  = flightAvg + (flightDist-flightAvg)/nstep;
  flightAvg2 = flightAvg2 + (flightDist*flightDist - flightAvg2)/nstep;

  // Compute change in distance every 10,000 steps
  G4double fltChange = flightAvg-lastFlight;
  lastFlight = flightAvg;

  // Compute change in RMS every 10,000 steps
  G4double RMS = sqrt(flightAvg2 - flightAvg*flightAvg);
  G4double RMSchange = RMS-lastRMS;
  lastRMS = RMS;

  // Scattering makes the path length longer, but only a factor of a few
  G4double pathScale = pathLen / flightDist;

  if (verboseLevel>1) {
    G4cout << " after " << nstep << " steps @ " << pos << G4endl;

    if (nstep%stepWindow == 1) G4cout << " new stepWindow" << G4endl;
    else {
      G4cout << " pos changed " << posShift << " since step "
	     << 1+((nstep-1)/stepWindow)*stepWindow << G4endl;
    }

    G4cout << " path " << pathLen << "  flight " << flightDist
	   << " : ratio " << pathScale << (pathScale>maxPathScale?" > ":" < ")
	   << maxPathScale << ")" << G4endl
	   << " flightAvg " << flightAvg << " changed by " << fltChange
	   << " RMS " << RMS << " changed by " << RMSchange << G4endl;
  }

  // Possible "stuck" conditions
  G4bool tooManySteps = (maxSteps>0 && nstep>maxSteps);
  G4bool excessPath = (pathScale > maxPathScale);
  G4bool windowFull = (nstep>=stepWindow && nstep%stepWindow == 0);
  G4bool samePos = (windowFull && posShift < minPosShift);
  G4bool notMoving = (windowFull && fabs(RMSchange) < minFlightRMS);

  if (verboseLevel>2) {
    G4cout << " tooManySteps " << tooManySteps;
    if (windowFull)
      G4cout << " samePos " << samePos << " notMoving " << notMoving;
    G4cout << G4endl;
  }

  return (tooManySteps || excessPath || samePos || notMoving);
}

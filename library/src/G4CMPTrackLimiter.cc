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

#include "G4CMPTrackLimiter.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPGeometryUtils.hh"
#include "G4CMPUtils.hh"
#include "G4ForceCondition.hh"
#include "G4ParticleChange.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VSolid.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include <limits.h>
#include <sstream>


// Only applies to G4CMP particles

G4bool G4CMPTrackLimiter::IsApplicable(const G4ParticleDefinition& pd) {
  return (G4CMP::IsPhonon(pd) || G4CMP::IsChargeCarrier(pd));
}


// Force killing if below cut

G4double G4CMPTrackLimiter::GetMeanFreePath(const G4Track& aTrack, G4double,
					    G4ForceCondition* condition) {
  G4bool changedLattice = UpdateMeanFreePathForLatticeChangeover(aTrack);
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

  return &aParticleChange;
}


// Evaluate current trackx

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
  G4double escape =
    ((postS->GetStepStatus() != fGeomBoundary) &&
     (postPV != GetCurrentVolume() || prePV != GetCurrentVolume()));

  if (verboseLevel>1) G4cout << " escape? " << escape << G4endl;
  
  return escape;
}

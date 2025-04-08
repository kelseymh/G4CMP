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


// Only applies to G4CMP particles

G4bool G4CMPTrackLimiter::IsApplicable(const G4ParticleDefinition& pd) {
  return (G4CMP::IsPhonon(pd) || G4CMP::IsChargeCarrier(pd));
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

  G4cout << std::setprecision(std::numeric_limits<double>::max_digits10) << G4endl;

  G4cout << "preStep Pos = " << step.GetPreStepPoint()->GetPosition() << G4endl;
  G4cout << "preStep status = " << step.GetPreStepPoint()->GetStepStatus() << G4endl;
  G4cout << "postStep Pos = " << step.GetPostStepPoint()->GetPosition() << G4endl;
  G4cout << "postStep status = " << step.GetPostStepPoint()->GetStepStatus() << G4endl;
  G4cout << "step dir = " << (step.GetPostStepPoint()->GetPosition() - step.GetPreStepPoint()->GetPosition()).unit() << G4endl;
  G4cout << "prePV Name = " << prePV->GetName() << G4endl;
  G4cout << "postPV Name = " << postPV->GetName() << G4endl;
  G4cout << "currPV Name = " << GetCurrentVolume()->GetName() << G4endl;

  G4VSolid* solid = GetCurrentVolume()->GetLogicalVolume()->GetSolid();
  EInside isIn = solid->Inside(step.GetPostStepPoint()->GetPosition());
  G4cout << "Value for surface: " << (isIn==kInside ? "inside" : isIn==kOutside ? "outside" : "surface") << G4endl;

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

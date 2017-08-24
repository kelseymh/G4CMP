/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPEnergyLimiter.cc
/// \brief Implementation of the G4CMPEnergyLimiter process, to kill tracks
///        falling below a minimum energy (set in G4CMPConfigManager).
//
// $Id$
//

#include "G4CMPEnergyLimiter.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPUtils.hh"
#include "G4Track.hh"
#include "G4ParticleChange.hh"
#include "G4ProcessType.hh"
#include <limits.h>


// Only applies to G4CMP particles

G4bool G4CMPEnergyLimiter::IsApplicable(const G4ParticleDefinition& pd) {
  return (G4CMP::IsPhonon(pd) || G4CMP::IsChargeCarrier(pd));
}


// Evaluate current track

G4bool G4CMPEnergyLimiter::BelowEnergyCut(const G4Track& track) const {
  G4double ecut =
    (G4CMP::IsChargeCarrier(track) ? G4CMPConfigManager::GetMinChargeEnergy()
     : G4CMP::IsPhonon(track) ? G4CMPConfigManager::GetMinPhononEnergy() : -1.);

  return (track.GetKineticEnergy() < ecut);
}


// Force killing if below cut

G4double G4CMPEnergyLimiter::
PostStepGetPhysicalInteractionLength(const G4Track& track, G4double,
                                     G4ForceCondition* condition) {
  *condition = (BelowEnergyCut(track) ? Forced : NotForced);
  return DBL_MAX;
}

G4VParticleChange* G4CMPEnergyLimiter::PostStepDoIt(const G4Track& track,
                                                    const G4Step& /*step*/) {
  aParticleChange.Initialize(track);

  if (BelowEnergyCut(track)) {
    aParticleChange.ProposeNonIonizingEnergyDeposit(track.GetKineticEnergy());
    aParticleChange.ProposeTrackStatus(fStopAndKill);
  }

  ClearNumberOfInteractionLengthLeft();		// All processes should do this!
  return &aParticleChange;
}

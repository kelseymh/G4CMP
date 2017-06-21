/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPTrackUtils.cc
/// \brief Free standing functions that perform operations on G4Tracks
///
//
// $Id$
//
// 201706211 M. Kelsey -- Non-templated utility functions

#include "G4CMPTrackUtils.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftTrackInfo.hh"
#include "G4CMPPhononTrackInfo.hh"
#include "G4CMPUtils.hh"
#include "G4CMPVTrackInfo.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4RandomDirection.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include "Randomize.hh"


// Use assigned particle type in track to set up initial TrackInfo

void G4CMP::AttachTrackInfo(const G4Track* track) {
  if (track) AttachTrackInfo(*track);
}

void G4CMP::AttachTrackInfo(const G4Track& track) {
  if (HasTrackInfo(track)) return;		// Don't replace existing!

  G4VPhysicalVolume* vol = track.GetVolume();
  G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();
  G4LatticePhysical* lat = LM->GetLattice(vol);
  
  if (IsPhonon(track)) {	// Assign random wavevector -- ?really right?
    auto trackInfo = new G4CMPPhononTrackInfo(lat, G4RandomDirection());
    AttachTrackInfo(track, trackInfo);
  }

  if (IsChargeCarrier(track)) {	// Assign random valley for electrons
    G4int valley = IsElectron(track) ? ChooseValley(lat) : -1;
    auto trackInfo = new G4CMPDriftTrackInfo(lat, valley);
    AttachTrackInfo(track, trackInfo);
  }
}


// Interrogate track to see if info already assigned (without type checking)

G4bool G4CMP::HasTrackInfo(const G4Track* track) {
  return (track && HasTrackInfo(*track));
}

G4bool G4CMP::HasTrackInfo(const G4Track& track) {
  G4VAuxiliaryTrackInformation* info = 
    track.GetAuxiliaryTrackInformation(G4CMPConfigManager::GetPhysicsModelID());

  return (nullptr != dynamic_cast<G4CMPVTrackInfo*>(info));
}

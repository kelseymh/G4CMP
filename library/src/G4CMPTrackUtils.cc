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
// 20170621 M. Kelsey -- Non-templated utility functions

#include "G4CMPTrackUtils.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftTrackInfo.hh"
#include "G4CMPGeometryUtils.hh"
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

  if (IsPhonon(track)) {
    AttachTrackInfo(track, G4RandomDirection());
  } else if (IsChargeCarrier(track)) {
    G4int valley = IsElectron(track) ? ChooseValley(GetLattice(track)) : -1;
    AttachTrackInfo(track, valley);
  }
}


// Create and initialize kinematics container for phonon track

void G4CMP::AttachTrackInfo(const G4Track* track, G4int valley) {
  if (track) AttachTrackInfo(*track, valley);
}

void G4CMP::AttachTrackInfo(const G4Track& track, G4int valley) {
  if (!IsChargeCarrier(track)) {
    G4Exception("G4CMP::AttachTrackInfo", "Utils003", JustWarning,
		"Cannot associate valley with phonon track.");
    AttachTrackInfo(track);		// Do default phonon attachment
    return;
  }

  AttachTrackInfo(track, new G4CMPDriftTrackInfo(GetLattice(track), valley));
}


// Create and initialize kinematics container for charged track

void G4CMP::AttachTrackInfo(const G4Track* track, const G4ThreeVector& kdir) {
  if (track) AttachTrackInfo(*track, kdir);
}

void G4CMP::AttachTrackInfo(const G4Track& track, const G4ThreeVector& kdir) {
  if (!IsPhonon(track)) {
    G4Exception("G4CMP::AttachTrackInfo", "Utils004", JustWarning,
		"Cannot associate wavevector with charged track.");
    AttachTrackInfo(track);		// Do default phonon attachment
    return;
  }
  
  AttachTrackInfo(track, new G4CMPPhononTrackInfo(GetLattice(track), kdir));
}


// Attach already created container to track, must be of G4CMP type

void G4CMP::AttachTrackInfo(const G4Track* track, G4CMPVTrackInfo* trackInfo) {
  if (track) AttachTrackInfo(*track, trackInfo);
}

void G4CMP::AttachTrackInfo(const G4Track& track, G4CMPVTrackInfo* trackInfo) {
  if (nullptr == trackInfo) return;

  track.SetAuxiliaryTrackInformation(G4CMPConfigManager::GetPhysicsModelID(),
				     trackInfo);
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


// Get physical lattice associated with track

G4LatticePhysical* G4CMP::GetLattice(const G4Track& track) {
  G4VPhysicalVolume* trkvol = track.GetVolume();
  if (!trkvol) trkvol = G4CMP::GetVolumeAtPoint(track.GetPosition());

  return G4LatticeManager::GetLatticeManager()->GetLattice(trkvol);
}

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
// 20190906 M. Kelsey -- Add function to look up process for track
// 20200829 M. Kelsey -- Don't override initial direction of phonons
// 20220907 G4CMP-316 -- Try using pre-step point to find lattice volume

#include "G4CMPTrackUtils.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftTrackInfo.hh"
#include "G4CMPGeometryUtils.hh"
#include "G4CMPPhononTrackInfo.hh"
#include "G4CMPUtils.hh"
#include "G4CMPVTrackInfo.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"

// Use assigned particle type in track to set up initial TrackInfo

void G4CMP::AttachTrackInfo(const G4Track* track) {
  if (track) AttachTrackInfo(*track);
}

void G4CMP::AttachTrackInfo(const G4Track& track) {
  if (HasTrackInfo(track)) return;		// Don't replace existing!

  if (IsPhonon(track)) {
    // Use track's initial momentum as wavevector direction
    AttachTrackInfo(track, track.GetMomentumDirection());
  } else if (IsChargeCarrier(track)) {
    G4int valley = IsElectron(track) ? ChooseValley(GetLattice(track)) : -1;
    AttachTrackInfo(track, valley);
  }else if (IsBogoliubovQP(track)){
      AttachTrackInfo(track, new G4CMPVTrackInfo(GetLattice(track)));
  }
}


// Create and initialize kinematics container for charged track

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


// Create and initialize kinematics container for phonon track

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
  //  G4cout << "REL-- back here, beginning of getlattice in trackutils." << G4endl;
  G4VPhysicalVolume* trkvol = track.GetVolume();
  if (!trkvol){
    trkvol = G4CMP::GetVolumeAtPoint(track.GetPosition());
    //    G4cout << "REL-- Have to use builtin G4CMP utils call to get volume at point... URK" << G4endl;
  }
  //G4cout << "REL-- physical volume at track point: " << trkvol->GetName() << G4endl;
  
  if (!G4LatticeManager::GetLatticeManager()->HasLattice(trkvol)
      && track.GetStep()) {
    trkvol = track.GetStep()->GetPreStepPoint()->GetPhysicalVolume();
  }

  return G4LatticeManager::GetLatticeManager()->GetLattice(trkvol);
}


// REL NB: I think this is no longer used -- keep for a little while to check but then delete!
// Need to define this for the case where we are transitioning from one lattice to another. Reason is that
// the above function uses track.GetVolume() whereas what needs to be called is track.GetNextVolume() for
// conditions where the lattices are not identical (which is I think exclusive to the needs of this tracked
// film response upgrade).
G4LatticePhysical* G4CMP::GetNextLattice(const G4Track & track)
{
  G4VPhysicalVolume* trkvol = track.GetNextVolume();
  if (!trkvol) trkvol = G4CMP::GetVolumeAtPoint(track.GetPosition());
  if (!G4LatticeManager::GetLatticeManager()->HasLattice(trkvol)
      && track.GetStep()) {
    //    G4cout << "REL running something not nice in TrackUtils: GetNextLattice." << G4endl;
    trkvol = track.GetStep()->GetPreStepPoint()->GetPhysicalVolume();
  }
  return G4LatticeManager::GetLatticeManager()->GetLattice(trkvol);
}


// Look up process by name associated with track

G4VProcess* G4CMP::FindProcess(const G4Track* track, const G4String& pname) {
  return (track ? G4CMP::FindProcess(*track, pname) : 0);
}

G4VProcess* G4CMP::FindProcess(const G4Track& track, const G4String& pname) {
  return G4CMP::FindProcess(track.GetDefinition(), pname);
}


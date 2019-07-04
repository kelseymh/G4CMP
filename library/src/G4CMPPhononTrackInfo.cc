/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPPhononTrackInfo.cc
/// \brief Implementation of the G4CMPPhononTrackInfo class. Used to store
/// auxiliary information that a G4Track can't store, but is necessary for
/// physics processes to know.
///
//
// $Id$
//
// 20161111 Initial commit - R. Agnese
// 20170728 M. Kelsey -- Replace "k" function args with "theK" (-Wshadow)

#include "G4CMPPhononTrackInfo.hh"

//G4Allocator<G4CMPPhononTrackInfo> G4CMPPhononTrackInfoAllocator;

G4CMPPhononTrackInfo::G4CMPPhononTrackInfo(const G4LatticePhysical* lat,
                                           G4ThreeVector theK)
  : G4CMPVTrackInfo(lat), waveVec(theK) {;}

void G4CMPPhononTrackInfo::Print() const {
//TODO
}

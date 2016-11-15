/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPDriftTrackInfo.cc
/// \brief Implementation of the G4CMPDriftTrackInfo class. Used to store
/// auxiliary information that a G4Track can't store, but is necessary for
/// physics processes to know.
///
//
// $Id$
//
// 20161111 Initial commit - R. Agnese

#include "G4CMPDriftTrackInfo.hh"
#include "G4LatticePhysical.hh"
#include "G4ParticleDefinition.hh"

//G4Allocator<G4CMPDriftTrackInfo> G4CMPDriftTrackInfoAllocator;

G4CMPDriftTrackInfo::G4CMPDriftTrackInfo(const G4LatticePhysical* lat,
                                         G4int valIdx) :
                                         G4CMPVTrackInfo(lat) {
  SetValleyIndex(valIdx);
}

void G4CMPDriftTrackInfo::SetValleyIndex(G4int valIdx) {
  // -1 is a valid value, so be careful coparing signed/unsigned.
  if (valIdx < -1 ||
      valIdx > static_cast<G4int>(Lattice()->NumberOfValleys() - 1)) {
    G4Exception("G4CMPDriftTrackInfo: Constructor", "DrfitTrackInfo001",
                EventMustBeAborted,
                "valley index parameter is out of bounds for the given lattice");
  }

  valleyIdx = valIdx;
}

void G4CMPDriftTrackInfo::Print() const {
//TODO
}

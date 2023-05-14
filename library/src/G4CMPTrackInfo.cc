/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPTrackInfo.cc
/// \brief Implementation of the G4CMPTrackInfo class. Used as a base class
/// for CMP particles to store auxiliary information that a G4Track can't
/// store, but is necessary for physics processes to know.
///
//
// $Id$
//
// 20161111 Initial commit - R. Agnese
// 20230514 M. Kelsey -- Rename G4CMPVTrackInfo to G4CMPTrackInfo (not virtual)

#include "G4CMPTrackInfo.hh"

G4CMPTrackInfo::G4CMPTrackInfo(const G4LatticePhysical* lat) :
  G4VAuxiliaryTrackInformation(), lattice(lat) {}

void G4CMPTrackInfo::Print() const {
//TODO
}

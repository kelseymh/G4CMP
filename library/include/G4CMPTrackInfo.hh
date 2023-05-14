/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPTrackInfo.hh
/// \brief Definition of the G4CMPTrackInfo class. Used as a base class
/// for CMP particles to store auxiliary information that a G4Track can't
/// store, but is necessary for physics processes to know.
///
//
// $Id$
//
// 20161111 Initial commit - R. Agnese
// 20230514 M. Kelsey -- Rename G4CMPVTrackInfo to G4CMPTrackInfo (not virtual)

#ifndef G4CMPTrackInfo_hh
#define G4CMPTrackInfo_hh 1

#include "G4VAuxiliaryTrackInformation.hh"

class G4LatticePhysical;

class G4CMPTrackInfo: public G4VAuxiliaryTrackInformation {
public:
  G4CMPTrackInfo() = delete;
  G4CMPTrackInfo(const G4LatticePhysical* lat);

  size_t ReflectionCount() const                           { return reflCount; }
  void IncrementReflectionCount()                               { ++reflCount; }

  const G4LatticePhysical* Lattice() const                   { return lattice; }
  void SetLattice(const G4LatticePhysical* lat)               { lattice = lat; }

  virtual void Print() const override;

private:
  size_t reflCount = 0; // Number of times track has been reflected
  const G4LatticePhysical* lattice; // The lattice the track is currently in
};

#endif

/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#include "G4CMPTrackInformation.hh"
#include "G4ThreeVector.hh"

// Default values are non-physical
G4CMPTrackInformation::G4CMPTrackInformation() :
  l0(0.), mc(0.), phononKVec(G4ThreeVector()), chargeValleyIdx(-1), reflCount(0)
{}

G4CMPTrackInformation::G4CMPTrackInformation(G4ThreeVector k) :
  l0(0.), mc(0.), phononKVec(k), chargeValleyIdx(-1), reflCount(0)
{}

G4CMPTrackInformation::G4CMPTrackInformation(G4double l,
                                             G4double m, G4int v) :
  l0(l), mc(m), phononKVec(G4ThreeVector()), chargeValleyIdx(v), reflCount(0)
{}

void G4CMPTrackInformation::Print() const
{
//TODO
}

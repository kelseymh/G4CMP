/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#include "G4CMPPhononTrackInfo.hh"

G4CMPPhononTrackInfo::G4CMPPhononTrackInfo(const G4LatticePhysical* lat,
                                           G4ThreeVector k) :
                                           G4CMPVTrackInfo(lat), waveVec(k) {}

void G4CMPPhononTrackInfo::Print() const {
//TODO
}

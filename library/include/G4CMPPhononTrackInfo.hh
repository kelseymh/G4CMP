/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPPhononTrackInfo.hh
/// \brief Definition of the G4CMPPhononTrackInfo class. Used to store
/// auxiliary information that a G4Track can't store, but is necessary for
/// physics processes to know.
///
//
// $Id$
//
// 20161111 Initial commit - R. Agnese
// 20170728 M. Kelsey -- Replace "k" function args with "theK" (-Wshadow)

#ifndef G4CMPPhononTrackInfo_hh
#define G4CMPPhononTrackInfo_hh 1

#include "G4CMPVTrackInfo.hh"
#include "G4ThreeVector.hh"

/* NOTE: Avoiding use of G4Allocator for performance reasons
#include "G4Allocator.hh"
class G4CMPPhononTrackInfo;
extern G4Allocator<G4CMPPhononTrackInfo> G4CMPPhononTrackInfoAllocator;
*/


class G4CMPPhononTrackInfo : public G4CMPVTrackInfo {
public:
  G4CMPPhononTrackInfo() = delete;
  G4CMPPhononTrackInfo(const G4LatticePhysical* lat, G4ThreeVector k);

/* NOTE: Avoiding use of G4Allocator for performance reasons
  void *operator new(size_t) noexcept {
    return static_cast<void*>(G4CMPPhononTrackInfoAllocator.MallocSingle());
  }

  void operator delete(void* info) noexcept {
    G4CMPPhononTrackInfoAllocator.FreeSingle(static_cast<G4CMPPhononTrackInfo*>(info));
  }
*/

  void SetK(G4ThreeVector theK)          { waveVec = theK; }
  void SetWaveVector(G4ThreeVector theK) { waveVec = theK; }
  G4ThreeVector k() const                { return waveVec; }
  G4ThreeVector WaveVector() const       { return waveVec; }

  virtual void Print() const override;

private:
  G4ThreeVector waveVec;
};

#endif

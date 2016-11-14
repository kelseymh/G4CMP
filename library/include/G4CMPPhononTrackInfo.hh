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

#ifndef G4CMPPhononTrackInfo_hh
#define G4CMPPhononTrackInfo_hh 1

#include "G4CMPVTrackInfo.hh"
#include "G4ThreeVector.hh"

class G4CMPPhononTrackInfo: public G4CMPVTrackInfo {
public:
  G4CMPPhononTrackInfo() = delete;
  G4CMPPhononTrackInfo(const G4LatticePhysical* lat, G4ThreeVector k);

  void SetK(G4ThreeVector k)                                    { waveVec = k; }
  void SetWaveVector(G4ThreeVector k)                           { waveVec = k; }
  G4ThreeVector k() const                                    { return waveVec; }
  G4ThreeVector WaveVector() const                           { return waveVec; }

  virtual void Print() const override;

private:
  G4ThreeVector waveVec;
};

#endif

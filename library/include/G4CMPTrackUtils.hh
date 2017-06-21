/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPTrackUtils.hh
/// \brief Free standing functions that perform operations on G4Tracks
///
//
// $Id$
//
// 20161111 Initial commit - R. Agnese
// 20170621 M. Kelsey -- Add non-templated utility functions

#include "globals.hh"

class G4CMPDriftTrackInfo;
class G4CMPPhononTrackInfo;
class G4CMPVTrackInfo;
class G4Track;

namespace G4CMP {
  void AttachTrackInfo(const G4Track* track);
  void AttachTrackInfo(const G4Track& track);

  template<class T> void AttachTrackInfo(const G4Track& track, T* trackInfo);

  template<class T> T* GetTrackInfo(const G4Track& track);

  G4bool HasTrackInfo(const G4Track* track);
  G4bool HasTrackInfo(const G4Track& track);
}

#include "G4CMPTrackUtils.icc"

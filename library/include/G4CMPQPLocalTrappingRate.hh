/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPQPLocalTrappingRate.hh
/// \brief Compute rate for QP local trapping
///
/// This is the rate class for the QP local trapping process. The
/// rate is drawn from lattice information rather than hardcoded.
//
// 20250922  G4CMP-219 -- First addition to this history (done at time
//                        of merge to develop)

#ifndef G4CMPQPLocalTrappingRate_hh
#define G4CMPQPLocalTrappingRate_hh 1

#include "G4CMPVScatteringRate.hh"


class G4CMPQPLocalTrappingRate : public G4CMPVScatteringRate {
public:
  G4CMPQPLocalTrappingRate() : G4CMPVScatteringRate("qpLocalTrapping") {;}
  virtual ~G4CMPQPLocalTrappingRate() {;}

  virtual G4double Rate(const G4Track& aTrack) const;
};

#endif	/* G4CMPQPLocalTrappingRate_hh */

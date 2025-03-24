/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPBogoliubovQPLocalTrappingRate.hh
/// \brief Compute rate for QP local trapping 
//
//


#ifndef G4CMPBogoliubovQPLocalTrappingRate_hh
#define G4CMPBogoliubovQPLocalTrappingRate_hh 1

#include "G4CMPVScatteringRate.hh"


class G4CMPBogoliubovQPLocalTrappingRate : public G4CMPVScatteringRate {
public:
  G4CMPBogoliubovQPLocalTrappingRate() : G4CMPVScatteringRate("BogoliubovQPLocalTrapping") {;}
  virtual ~G4CMPBogoliubovQPLocalTrappingRate() {;}

  virtual G4double Rate(const G4Track& aTrack) const;
};

#endif	/* G4CMPBogoliubovQPLocalTrappingRate_hh */

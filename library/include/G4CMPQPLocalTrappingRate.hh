/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPQPLocalTrappingRate.hh
/// \brief Compute rate for QP local trapping 
//
//


#ifndef G4CMPQPLocalTrappingRate_hh
#define G4CMPQPLocalTrappingRate_hh 1

#include "G4CMPVScatteringRate.hh"


class G4CMPQPLocalTrappingRate : public G4CMPVScatteringRate {
public:
  G4CMPQPLocalTrappingRate() : G4CMPVScatteringRate("QPLocalTrapping") {;}
  virtual ~G4CMPQPLocalTrappingRate() {;}

  virtual G4double Rate(const G4Track& aTrack) const;
};

#endif	/* G4CMPQPLocalTrappingRate_hh */

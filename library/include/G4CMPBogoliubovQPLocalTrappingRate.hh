/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPPhononScatteringRate.hh
/// \brief Compute rate for phonon impurity scattering (mode mixing)
//
// $Id$
//
// 20170815  Move G4CMPProcessUtils inheritance to base class

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

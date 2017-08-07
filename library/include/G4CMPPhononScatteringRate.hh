/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPPhononScatteringRate.hh
/// \brief Compute rate for phonon impurity scattering (mode mixing)
//
// $Id$

#ifndef G4CMPPhononScatteringRate_hh
#define G4CMPPhononScatteringRate_hh 1

#include "G4CMPVScatteringRate.hh"
#include "G4CMPProcessUtils.hh"


class G4CMPPhononScatteringRate : public G4CMPVScatteringRate,
				  public G4CMPProcessUtils {
public:
  G4CMPPhononScatteringRate() : G4CMPVScatteringRate("PhononScattering") {;}
  virtual ~G4CMPPhononScatteringRate() {;}

  virtual G4double Rate(const G4Track& aTrack) const;
};

#endif	/* G4CMPPhononScatteringRate_hh */

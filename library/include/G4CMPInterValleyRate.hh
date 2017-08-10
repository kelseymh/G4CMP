/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPInterValleyRate.hh
/// \brief Compute electron intervalley scattering rate using matrix elements
//
// $Id$

#ifndef G4CMPInterValleyRate_hh
#define G4CMPInterValleyRate_hh 1

#include "G4CMPVScatteringRate.hh"
#include "G4CMPProcessUtils.hh"


class G4CMPInterValleyRate : public G4CMPVScatteringRate,
			     public G4CMPProcessUtils {
public:
  G4CMPInterValleyRate() : G4CMPVScatteringRate("InterValley") {;}
  virtual ~G4CMPInterValleyRate() {;}

  virtual G4double Rate(const G4Track& aTrack) const;
};

#endif	/* G4CMPInterValleyRate_hh */

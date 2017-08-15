/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPDownconversionRate.hh
/// \brief Compute rate for phonon anharmonic decay (downconversion)
//
// $Id$
//
// 20170815  Move G4CMPProcessUtils inheritance to base class

#ifndef G4CMPDownconversionRate_hh
#define G4CMPDownconversionRate_hh 1

#include "G4CMPVScatteringRate.hh"


class G4CMPDownconversionRate : public G4CMPVScatteringRate {
public:
  G4CMPDownconversionRate() : G4CMPVScatteringRate("Downconversion") {;}
  virtual ~G4CMPDownconversionRate() {;}

  virtual G4double Rate(const G4Track& aTrack) const;
};

#endif	/* G4CMPDownconversionRate_hh */

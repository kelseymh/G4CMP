/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPLukeEmissionRate.hh
/// \brief Compute emission rate for Luke-Neganov phonons.
//
// $Id$
//
// 20170815  Move G4CMPProcessUtils inheritance to base class
// 20170907  Make process non-forced; TimeStepper will trigger recalculation

#ifndef G4CMPLukeEmissionRate_hh
#define G4CMPLukeEmissionRate_hh 1

#include "G4CMPVScatteringRate.hh"


class G4CMPLukeEmissionRate : public G4CMPVScatteringRate {
public:
  G4CMPLukeEmissionRate() : G4CMPVScatteringRate("Luke") {;}
  virtual ~G4CMPLukeEmissionRate() {;}

  virtual G4double Rate(const G4Track& aTrack) const;
};

#endif	/* G4CMPLukeEmissionRate_hh */

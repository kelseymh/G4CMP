/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPIVRateLinear.hh
/// \brief Compute electron intervalley scattering rate using Linear
///	   parametrization vs. electric field.
//
// $Id$
//
// 20170815  Move G4CMPProcessUtils inheritance to base class

#ifndef G4CMPIVRateLinear_hh
#define G4CMPIVRateLinear_hh 1

#include "G4CMPVScatteringRate.hh"


class G4CMPIVRateLinear : public G4CMPVScatteringRate {
public:
  G4CMPIVRateLinear() : G4CMPVScatteringRate("IVLinear") {;}
  virtual ~G4CMPIVRateLinear() {;}

  virtual G4double Rate(const G4Track& aTrack) const;
};

#endif	/* G4CMPIVRateLinear_hh */

/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPIVRateQuadratic.hh
/// \brief Compute electron intervalley scattering rate using Quadratic
///	   parametrization vs. electric field.
//
// $Id$
//
// 20170815  Move G4CMPProcessUtils inheritance to base class
// 20240823  Change name to plain "Quadratic", to match UseRateModel()
// 20250515  Apply IV energy thresholds to reduce zero-voltage scatters

#ifndef G4CMPIVRateQuadratic_hh
#define G4CMPIVRateQuadratic_hh 1

#include "G4CMPVScatteringRate.hh"


class G4CMPIVRateQuadratic : public G4CMPVScatteringRate {
public:
  G4CMPIVRateQuadratic() : G4CMPVScatteringRate("Quadratic") {;}
  virtual ~G4CMPIVRateQuadratic() {;}

  virtual G4double Rate(const G4Track& aTrack) const;
  virtual G4double Threshold(G4double Eabove=0.) const;
};

#endif	/* G4CMPIVRateQuadratic_hh */

/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPKaplanQPRate.hh
/// \brief Compute Rates for QP scattering and Recombination
//
//
//

#ifndef G4CMPQPScatteringRate_hh
#define G4CMPQPScatteringRate_hh 1

#include "G4CMPVScatteringRate.hh"

class G4CMPQPScatteringRate : public G4CMPVScattering {
public:
  G4CMPQPScatteringRate() :  G4CMPVScatteringRate("Bogoliubov") {;}

  virtual ~ G4CMPQPScatteringRate() {;}

  virtual G4double Rate(const G4Track& aTrack) const;

  virtual Threshold(G4double energy=0.) const;

  virtual void LoadDataForTrack(const G4Track* track);

protected:

private:

};

#endif

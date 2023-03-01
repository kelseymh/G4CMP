/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// \file library/include/G4CMPQPRecombinationRate.hh
// \brief Compute rate for quasiparticle scattering


#ifndef G4CMPQPRecombinationRate_hh
#define G4CMPQPRecombinationRate_hh 1

#include "G4CMPVScatteringRate.hh"

class G4CMPQPRecombinationRate : public G4CMPVScatteringRate {
public:
  G4CMPQPRecombinationRate() : G4CMPVScatteringRate("QPRecombination") {;}
  virtual ~G4CMPQPRecombinationRate() {;}


  virtual G4double Rate(const G4Track& aTrack) const;
  virtual G4double Rate(G4double energy, G4double temp=0) const;
private:
  //helper functions
  G4double Recombination(G4double qp_E, G4double x, G4double T, G4double Tc, G4double Gap0, G4double b);
  G4double Numeric_Integral(N_emm G4double qp_E,  G4double T, G4double Tc, G4double Gap0, G4double b);

}

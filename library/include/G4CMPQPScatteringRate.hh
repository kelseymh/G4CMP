/***********************************************************************\
 *  This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/
// \file library/include/G4CMPQPScatteringRate.hh
// \brief Compute rate for quasiparticle scattering


#ifndef G4CMPQPScatteringRate_hh
#define G4CMPQPScatteringRate_hh 1


#include "G4CMPVScatteringRate.hh"

class G4CMPQPScatteringRate : public G4CMPVScatteringRate {
public:
  G4CMPQPScatteringRate() : G4CMPVScatteringRate("QPScattering") {;}
  virtual ~G4CMPQPScatteringRate() {;}


  virtual G4double Rate(const G4Track& aTrack) const;
  virtual G4double Rate(G4double energy, G4double temp=0) const;
private:

  //helper functions
  //emmision integrand
  G4double Scattering_Emmision(G4double qp_E, G4double x, G4double T, G4double Tc, G4double Gap0, G4double b);
  //absorptin integrand
  G4double Scattering_Absoption(G4double qp_E, G4double x, G4double T, G4double Tc, G4double Gap0, G4double b);
  G4double Numeric_Integral(int N_emm, int N_abs, G4double qp_E,  G4double T, G4double Tc, G4double Gap0, G4double b);
}


#endif /*G4CMPQPScatteringRate*/

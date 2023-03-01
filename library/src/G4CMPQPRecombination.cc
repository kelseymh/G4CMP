/***********************************************************************\
 *This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/
// \file library/src/G4CMPQPRecombinationRate.cc
// \brief Compute rate for quasiparticle recombination

// $Id$


#include "G4CMPQPScatteringRate.hh"
#include "G4PhysicalConstants.hh"
#include "G4Track.hh"
#include "G4CMPUtils.hh"
#include "G4CMPConfigManager.hh"
#include <cmath>

//Global variables for Al, these should be moves into the crystal maps
G4double Z = 1.43; //renormalization paramter Z(0)
G4double b = 0.000317; //meV ^-2
G4double Tc = 1.1810; //K
G4double Gap0 = .00017 // eV


//Recombination rate for tracks
G4double G4CMPQPScatteringRate::Rate(const G4Track& aTrack) const {

  G4double qp_E = GetKineticEnergy(aTrack);
  G4double temperature = G4CMPConfigManager::GetTemperature();
  G4double coeff = 2*pi/hbar_Planck / Z / (1-Fermi(temperature, qp_E));
   
  return coeff * Numeric_Integral(10000, 10000, qp_E, temperature,  Tc, Gap0, b);
} 

//Scattering Rate for energy and temp
G4double G4CMPQPScatteringRate::Rate(G4double temperature, G4double energy) const {

  G4double coeff = 2*pi/hbar_Planck / Z / (1-Fermi(temperature, energy));
  return coeff * Numeric_Integral(10000, energy, temperature,  Tc, Gap0, b);
}

//From kaplan 1976 eq 8
G4double G4CMPQPRecombinationRate::Recombination(G4double qp_E, G4double x, G4double T, G4double Tc, G4double Gap0, G4double b) { 
//qp energy, integration var, temperature of qp, critical temp,  chracateristic constant of material
  G4complex<G4double> comp_part(  ((x - qp_E)*(x - qp_E)) / ((x - qp_E)*(x - qp_E) - Gap(T,Gap0,Tc)*Gap(T,Gap0,Tc)), 0);
  G4double comp_real = sqrt(comp_part).real();

  return b*x*x* comp_real * (1 - Gap(T, Gap0, Tc)*Gap(T, Gap0, Tc)/qp_E/(x - qp_E)) * (Bose(T, x) + 1) *  Fermi(T, x - qp_E);
}

G4double G4CMPQPScatteringRate::Numeric_Integral(int N, G4double qp_E,  G4double T, G4double Tc, G4double Gap0, G4double b) {
  //Integral from a to b
  G4double a = qp_E + Gap(T, Gap0, Tc);
  G4double b = 1000; //This integral should converge quickly and if we are on order of 1000 eV, we have bigger problems

  G4double sum = Recombination(qp_E, a, T, Tc, Gap0, b)/2 + Recombination(qp_E, b, T, Tc, Gap0, b)/2;
  
  G4double var_x = 0;
  for(int i = 1; i < N-1; i++)  {
    var_x = a + i*(b-a)/N;
    sum_emm += Recombination(qp_E, var_x, T, Tc, Gap0, b);
  }

  sum*= (b-a)/N;

  return sum;
}

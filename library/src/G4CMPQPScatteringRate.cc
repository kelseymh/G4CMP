/// \file library/src/G4CMPQPScatteringRate.cc
///// \brief Compute rate for quasiparticle scattering




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

//Scattering rate for tracks
G4double G4CMPQPScatteringRate::Rate(const G4Track& aTrack) const {

  G4double qp_E = GetKineticEnergy(aTrack);
  G4double temperature = G4CMPConfigManager::GetTemperature();
  G4double coeff = 2*pi/hbar_Planck / Z / (1-Fermi(temperature, qp_E));
   
  return coeff * Numeric_Integral(10000, 10000, qp_E, temperature,  Tc, Gap0, b);
} 

//Scattering Rate for energy and temp
G4double G4CMPQPScatteringRate::Rate(G4double temperature, G4double energy) const {

  G4double coeff = 2*pi/hbar_Planck / Z / (1-Fermi(temperature, energy));
  return coeff * Numeric_Integral(10000, 10000, energy, temperature,  Tc, Gap0, b);
}

//From kaplan 1976 eq 8
G4double G4CMPQPScatteringRate::Scattering_Emmision(G4double qp_E, G4double x, G4double T, G4double Tc, G4double Gap0, G4double b) { 
//qp energy, integration var, temperature of qp, critical temp,  chracateristic constant of material
  G4complex<G4double> comp_part(  ((qp_E - x)*(qp_E - x)) / ((qp_E - x)*(qp_E - x) - Gap(T,Gap0,Tc)*Gap(T,Gap0,Tc)), 0);
  G4double comp_real = sqrt(comp_part).real(); 

 return b*x*x* comp_real * (1 - Gap(T, Gap0, Tc)*Gap(T, Gap0, Tc)/qp_E/(qp_E - x)) * (Bose(T, x) + 1) * (1 - Fermi(T, qp_E - x));
 
}
G4double Scattering_Absorption(G4double qp_E, G4double x, G4double T, G4double Tc, G4double Gap0, G4double b) {
  G4complex<G4double> comp_part(  ((qp_E + x)*(qp_E + x)) / ((qp_E + x)*(qp_E + x) - Gap(T,Gap0,Tc)*Gap(T,Gap0,Tc)), 0);
  G4double comp_real = sqrt(comp_part).real();

  return b*x*x* comp_real * (1 - Gap(T, Gap0, Tc)*Gap(T, Gap0, Tc)/qp_E/(qp_E + x)) * Bose(T, x) * (1 - Fermi(T, qp_E + x));
}

G4double G4CMPQPScatteringRate::Numeric_Integral(int N_emm, int N_abs, G4double qp_E,  G4double T, G4double Tc, G4double Gap0, G4double b) {
  //Integral from 0 to b
  G4double b_emm = qp_E - Gap(T, Gap0, Tc);
  G4double b_abs = 1000.; //if we have energies at 1000 ev, we have bigger problems, and this integral should quickly converge

  G4double sum_abs = Scattering_Absorption(qp_E, b_abs, T, Tc, Gap0, b)/2; //sum for absoption
  G4double sum_emm = Scattering_Emmision(qp_E, b_emm, T, Tc, Gap0, b)/2; //sum for emmision

  G4double var_x = 0;
  for(int i = 1; i < N_emm-1; i++)  {
    var_x = i*b_emm/N_emm;
    sum_emm += Scattering_Emmision(qp_E, var_x, T, Tc, Gap0, b);
  }    
  for(int i = 1; i < N_abs-1; i++)  {
    var_x = i*b_abs/N_abs;
    sum_emm += Scattering_Absorption(qp_E, var_x, T, Tc, Gap0, b);
  }

  sum_emm*= b_emm/N_emm;
  sum_abs*= b_abs/Nabs;

  return sum_emm+sum_abs;
}


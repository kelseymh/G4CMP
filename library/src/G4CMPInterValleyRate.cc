/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPInterValleyRate.cc
/// \brief Compute electron intervalley scattering rate using matrix elements
//
// $Id$

#include "G4CMPInterValleyRate.hh"
#include "G4Field.hh"
#include "G4FieldManager.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalVolume.hh"
#include "G4PhysicalConstants.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include <math.h>


// Scattering rate is computed from electric field

G4double G4CMPInterValleyRate::Rate(const G4Track& aTrack) const {
  const_cast<G4CMPInterValleyRate*>(this)->LoadDataForTrack(&aTrack);

  // Track kinematics
  G4double velocity = GetVelocity(aTrack);
  G4double energy = GetKineticEnergy(aTrack);

  // Crystal properties (will want to get from lattice, material, etc.)
  G4double rho =
    aTrack.GetVolume()->GetLogicalVolume()->GetMaterial()->GetDensity();

  G4double temp = 0.015*kelvin;		// Set in lattice? material?
  G4double kT = k_Boltzmann * temp;

  G4double hbar_sq  = hbar_Planck * hbar_Planck;
  G4double hbar_4th = hbar_sq * hbar_sq;

  G4double mass_electron = electron_mass_c2/c_squared;  
  G4double ml = theLattice->GetMassTensor().xx();	// m(parallel)
  G4double mt = theLattice->GetMassTensor().yy();	// m(perp)

  G4double M_D = cbrt(ml*mt*mt);    	// Average carrier mass
  G4double M_D3half = sqrt(M_D*M_D*M_D);
 
  //Acoustic Phonon Scattering 
  G4double D_ac  = 11*eV;            
  G4double alpha = 0.3/eV;                     

  // Useful constants for expressions below
  G4double D_ac_sq = D_ac*D_ac;
  G4double velocity_sq = velocity*velocity;

  G4double arate = ( sqrt(2)*kT * M_D3half * D_ac_sq *
		     sqrt(energy*(1+alpha*energy))*(1+2*alpha*energy)
		     / (pi*hbar_4th*rho*velocity_sq) );

  if (verboseLevel>2) G4cout << " Acoustic phonon rate " <<  arate << G4endl;

  // Optical phonon Scattering equation
  G4double D_op[]  = { 3e10*eV, 2e9*eV };	// Are these values correct?
  G4double hw_op[] = { 27.3e-3*eV, 10.3e-3*eV }; // Shouldn't w_op be frequency?
  G4double orate[]  = { 0,0 };
  G4double orateTotal = 0;
  
  for (int i = 0; i<2; i++) {
    if (energy <= hw_op[i]) continue;		// Apply threshold behaviour

    G4double energy_minus_hw_op = energy - hw_op[i];
    G4double D_op_sq = D_op[i]*D_op[i];
    G4double alpha_times_ehw_op = alpha * energy_minus_hw_op;
    G4double everything_under_sqrt =
      sqrt(energy_minus_hw_op*(1 + alpha_times_ehw_op));

    orate[i] = kT * M_D3half * D_op_sq * everything_under_sqrt *
      (1 + 2*alpha_times_ehw_op) / (sqrt(2) * pi* hbar_sq *rho*hw_op[i]);

    if (verboseLevel>2)
      G4cout << " Optical phonon rate [" << i << "] " << orate[i] << G4endl;

    orateTotal += orate[i];
  }

  if (verboseLevel>2) G4cout << " Optical Phonon rate " << orateTotal << G4endl;
 
  //Neutral Impurities
  G4double n_I = 1e11/cm3; 			// Number density of inpurities
  G4double epsilon_r = 16.2;			// Relative permittivity of Ge
  G4double E_T = (M_D/mass_electron) / epsilon_r;
  
  G4double nrate = (4*sqrt(2)* n_I * hbar_sq * sqrt(energy)) /
    ( M_D3half * (energy + E_T));

  if (verboseLevel>2) G4cout << " Neutral Impurities " << nrate << G4endl;

  G4double rate = nrate + orateTotal + arate;
  if (verboseLevel>1) G4cout << "IV rate = " << rate/hertz << " Hz" << G4endl;
  return rate;
}

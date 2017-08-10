/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPInterValleyRate.cc
/// \brief Compute electron intervalley scattering rate using matrix elements
//
// $Id$

#include "G4CMPInterValleyRate.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalVolume.hh"
#include "G4PhysicalConstants.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include <math.h>


// Scattering rate is computed from matrix elements

G4double G4CMPInterValleyRate::Rate(const G4Track& aTrack) const {
  const_cast<G4CMPInterValleyRate*>(this)->LoadDataForTrack(&aTrack);

  // Track kinematics
  G4double velocity = GetVelocity(aTrack);
  G4double energy = GetKineticEnergy(aTrack);

  // Crystal properties (will want to get from lattice, material, etc.)
  G4double density = theLattice->GetDensity();

  G4double temp = 0.015*kelvin;		// Set in lattice? material?
  G4double kT = k_Boltzmann * temp;

  G4double hbar_sq  = hbar_Planck * hbar_Planck;
  G4double hbar_4th = hbar_sq * hbar_sq;

  G4double m_electron = electron_mass_c2/c_squared;  
  G4double m_DOS = theLattice->GetElectronDOSMass();
  G4double m_DOS3half = sqrt(m_DOS*m_DOS*m_DOS);
 
  //Acoustic Phonon Scattering 
  G4double D_ac  = theLattice->GetAcousticDeform();
  G4double alpha = theLattice->GetAlpha();	// Non-parabolicity scale

  // Useful constants for expressions below
  G4double D_ac_sq = D_ac*D_ac;
  G4double velocity_sq = velocity*velocity;

  G4double arate = ( sqrt(2)*kT * m_DOS3half * D_ac_sq *
		     sqrt(energy*(1+alpha*energy))*(1+2*alpha*energy)
		     / (pi*hbar_4th*density*velocity_sq) );

  if (verboseLevel>2) G4cout << " Acoustic phonon rate " <<  arate << G4endl;

  // Optical Phonon Scattering
  G4double orateTotal = 0;

  G4int N_op = theLattice->GetNOptical();
  for (G4int i = 0; i<N_op; i++) {
    G4double D_op = theLattice->GetOpticalDeform(i);
    G4double Emin_op = theLattice->GetOpticalEnergy(i);

    if (energy <= Emin_op) continue;		// Apply threshold behaviour

    G4double dE = energy - Emin_op;		// Energy above threshold
    G4double D_op_sq = D_op*D_op;

    G4double orate = ( kT * m_DOS3half * D_op_sq *
		       sqrt(dE*(1 + alpha*dE))*(1 + 2*alpha*dE)
		       / (sqrt(2)*pi*hbar_sq*density*Emin_op) );

    if (verboseLevel>2)
      G4cout << " Optical phonon rate [" << i << "] " << orate << G4endl;

    orateTotal += orate;
  }

  if (verboseLevel>2) G4cout << " Optical Phonon rate " << orateTotal << G4endl;
 
  //Neutral Impurities
  G4double n_I = theLattice->GetImpurities();		// Number density
  G4double epsilon_r = theLattice->GetPermittivity();	// Relative
  G4double E_T = 0.75*eV * (m_DOS/m_electron) / epsilon_r;
  
  G4double nrate = ( 4*sqrt(2)* n_I * hbar_sq * sqrt(energy)
		     / (m_DOS3half * (energy+E_T)) );

  if (verboseLevel>2) G4cout << " Neutral Impurities " << nrate << G4endl;

  G4double rate = nrate + orateTotal + arate;
  if (verboseLevel>1) G4cout << "IV rate = " << rate/hertz << " Hz" << G4endl;
  return rate;
}

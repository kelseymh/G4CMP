/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPInterValleyRate.cc
/// \brief Compute electron intervalley scattering rate using matrix elements
//
// $Id$
//
// 20170821  Follow Aubry-Fortuna (2005) for separate D0 and D1 scattering

#include "G4CMPInterValleyRate.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalVolume.hh"
#include "G4PhysicalConstants.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include <math.h>


// Initialize lattice parameters used in matrix element calculations

void G4CMPInterValleyRate::LoadDataForTrack(const G4Track* track) {
  G4CMPProcessUtils::LoadDataForTrack(track);

  // Should temperature be a lattice configuration?
  kT = k_Boltzmann * 0.015*kelvin;

  density = theLattice->GetDensity();
  alpha = theLattice->GetAlpha();
  nValley = theLattice->NumberOfValleys()-1;

  m_DOS = theLattice->GetElectronDOSMass();
  m_DOS3half = sqrt(m_DOS*m_DOS*m_DOS);
}


// Scattering rate is computed from matrix elements

G4double G4CMPInterValleyRate::Rate(const G4Track& aTrack) const {
  const_cast<G4CMPInterValleyRate*>(this)->LoadDataForTrack(&aTrack);

  // Initialize numerical buffers
  vTrk = GetVelocity(aTrack);		  // Track kinematics
  eTrk = GetKineticEnergy(aTrack);

  G4double arate = acousticRate();
  if (verboseLevel>2) G4cout << " Acoustic phonon rate " <<  arate << G4endl;

  G4double orate = opticalD0Rate() + opticalD1Rate();
  if (verboseLevel>2) G4cout << " Optical Phonon rate " << orate << G4endl;
 
  G4double nrate = scatterRate();
  if (verboseLevel>2) G4cout << " Neutral Impurities " << nrate << G4endl;

  G4double rate = nrate + orate + arate;
  if (verboseLevel>1) G4cout << "IV rate = " << rate/hertz << " Hz" << G4endl;
  return rate;
}


// Compute components of overall intervalley rate

G4double G4CMPInterValleyRate::acousticRate() const {
  G4double D_ac  = theLattice->GetAcousticDeform();
  G4double D_ac_sq = D_ac*D_ac;
  G4double vTrk_sq = vTrk*vTrk;

  return ( sqrt(2)*kT * m_DOS3half * D_ac_sq *
	   sqrt(eTrk*(1+alpha*eTrk))*(1+2*alpha*eTrk)
	   / (pi*hbar_4th*density*vTrk_sq) );
}

G4double G4CMPInterValleyRate::opticalD0Rate() const {
  G4double totalD0 = 0.;

  G4int N_op = theLattice->GetNOptical(0);
  for (G4int i = 0; i<N_op; i++) {
    G4double Emin_op = theLattice->GetOpticalEnergy(0,i);
    if (eTrk <= Emin_op) continue;		// Apply threshold behaviour

    G4double dE = eTrk - Emin_op;		// Energy above threshold
    G4double D_op = theLattice->GetOpticalDeform(0,i);
    G4double D_op_sq = D_op*D_op;

    G4double orate = ( nValley * kT * m_DOS3half * D_op_sq *
		       sqrt(dE*(1 + alpha*dE))*(1 + 2*alpha*dE)
		       / (sqrt(2)*pi*hbar_sq*density*Emin_op) );

    if (verboseLevel>2)
      G4cout << " D0 phonon rate [" << i << "] " << orate << G4endl;

    totalD0 += orate;
  }

  return totalD0;
}

G4double G4CMPInterValleyRate::opticalD1Rate() const {
  G4double totalD1 = 0.;

  G4int N_op = theLattice->GetNOptical(1);
  for (G4int i = 0; i<N_op; i++) {
    G4double Emin_op = theLattice->GetOpticalEnergy(1,i);
    if (eTrk <= Emin_op) continue;		// Apply threshold behaviour

    G4double dE = eTrk - Emin_op;		// Energy above threshold
    G4double D_op = theLattice->GetOpticalDeform(1,i);
    G4double D_op_sq = D_op*D_op;

    G4double orate = ( sqrt(2)*nValley * m_DOS*m_DOS3half * D_op_sq *
		       sqrt(dE*(1 + alpha*dE))*(1 + 2*alpha*dE) *
		       (dE*(1+alpha*dE) + (1+alpha*eTrk))
		       / (pi*hbar_4th*density*Emin_op) );

    if (verboseLevel>2)
      G4cout << " D1 phonon rate [" << i << "] " << orate << G4endl;

    totalD1 += orate;
  }

  return totalD1;
}

G4double G4CMPInterValleyRate::scatterRate() const {
  G4double n_I = theLattice->GetImpurities();		// Number density
  G4double epsilon_r = theLattice->GetPermittivity();	// Relative
  G4double E_T = 0.75*eV * (m_DOS/m_electron) / epsilon_r;
  
  return ( 4*sqrt(2)* n_I * hbar_sq * sqrt(eTrk)
	   / (m_DOS3half * (eTrk+E_T)) );
}

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
// 20170830  Follow Jacoboni, with unified D0/D1 expression and units; drop
//		acoustic rate, as it is _intra_valley.
// 20170919  Add interface for threshold identification

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

  uSound = (2.*theLattice->GetTransverseSoundSpeed()
	    + theLattice->GetSoundSpeed()) / 3.;

  density = theLattice->GetDensity();
  alpha = theLattice->GetAlpha();
  nValley = 2*theLattice->NumberOfValleys()-1;		// From symmetry

  m_DOS = theLattice->GetElectronDOSMass();
  m_DOS3half = sqrt(m_DOS*m_DOS*m_DOS);
}


// Scattering rate is computed from matrix elements

G4double G4CMPInterValleyRate::Rate(const G4Track& aTrack) const {
  const_cast<G4CMPInterValleyRate*>(this)->LoadDataForTrack(&aTrack);

  // Initialize numerical buffers
  eTrk = GetKineticEnergy(aTrack);
  if (verboseLevel>1)
    G4cout << "G4CMPInterValleyRate eTrk " << eTrk/eV << " eV" << G4endl;

  G4double orate = opticalRate();
  if (verboseLevel>2) G4cout << "IV phonons  " << orate/hertz << " Hz" << G4endl;
 
  G4double nrate = scatterRate();
  if (verboseLevel>2) G4cout << "IV neutrals " << nrate/hertz << " Hz" << G4endl;

  G4double rate = nrate + orate;
  if (verboseLevel>1) G4cout << "IV rate = " << rate/hertz << " Hz" << G4endl;
  return rate;
}


// Compute components of overall intervalley rate

G4double G4CMPInterValleyRate::acousticRate() const {
  G4double D_ac  = theLattice->GetAcousticDeform();
  G4double D_ac_sq = D_ac*D_ac;

  return ( sqrt(2)*kT * m_DOS3half * D_ac_sq * energyFunc(eTrk)
	   / (pi*hbar_4th*density*uSound*uSound) );
}

G4double G4CMPInterValleyRate::opticalRate() const {
   // FIXME:  Rate should not have 'kT', but leaving it out ruins drift curve
  G4double scale = nValley*/*kT**/m_DOS3half / (sqrt(2)*pi*hbar_sq*density);

  G4double total = 0.;
  G4int N_op = theLattice->GetNIVDeform();
  for (G4int i = 0; i<N_op; i++) {
    G4double Emin_op = theLattice->GetIVEnergy(i);
    if (eTrk <= Emin_op) continue;		// Apply threshold behaviour

    G4double D_op = theLattice->GetIVDeform(i);
    G4double oscale = scale * D_op*D_op / Emin_op;

    G4double Efunc = energyFunc(eTrk-Emin_op);	// Energy above threshold

    G4double orate = oscale * Efunc;

    if (verboseLevel>2) {
      G4cout << " oscale[" << i << "] " << oscale << " Efunc " << Efunc
	     << "\n phonon rate [" << i << "] " << orate/hertz << " Hz"
	     << G4endl;
    }

    total += orate;
  }

  return total;
}

G4double G4CMPInterValleyRate::scatterRate() const {
  G4double n_I = theLattice->GetImpurities();		// Number density
  G4double epsilon_r = theLattice->GetPermittivity();	// Dielectric constant
  G4double E_T = 0.75*eV * (m_DOS/m_electron) / epsilon_r;
  
  return ( 4*sqrt(2)* n_I * hbar_sq * sqrt(eTrk)
	   / (m_DOS3half * (eTrk+E_T)) );
}


// Identify next energy threshold (if any) above specified input

G4double G4CMPInterValleyRate::Threshold(G4double Eabove) const {
  // Get list of all energy thresholds and sort
  std::vector<G4double> E_op = theLattice->GetIVEnergy();
  if (E_op.empty()) return 0.;

  // Put energies in order, find nearest entry above input value
  std::sort(E_op.begin(), E_op.end());
  std::vector<G4double>::const_iterator thresh =
    std::upper_bound(E_op.begin(), E_op.end(), Eabove);

  return (thresh == E_op.end() ? 0. : *thresh);
}

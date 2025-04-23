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
// 20250423  Add parabolic ellipsoidal bands IV scattering rate.

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

  density = theLattice->GetDensity();
  alpha = theLattice->GetAlpha();
  m_DOS = theLattice->GetElectronDOSMass();
  m_DOS3half = sqrt(m_DOS*m_DOS*m_DOS);

}


// Scattering rate from first principles

G4double G4CMPInterValleyRate::Rate(const G4Track& aTrack) const {
  const_cast<G4CMPInterValleyRate*>(this)->LoadDataForTrack(&aTrack);
    
  // Initialize Energy and momentum
  eTrk = GetKineticEnergy(aTrack);    
  ivalley = GetValleyIndex(aTrack);
  ptrk = GetLocalMomentum(aTrack);
  ktrk = theLattice->MapPtoK(ivalley, ptrk);
  kHV = theLattice->EllipsoidalToSphericalTranformation(ivalley, ktrk);
  kmag = kHV.mag();  
  
  IVprob = {};		// Store IV rates
  G4double totalIVRate = 0.;      // Total IV rate
  G4int N_op = theLattice->GetNIVDeform();		// # of IV transitions possible

  // Going through each phonon mode
  for (G4int i = 0; i<N_op; i++) {
      
    G4double Emin_iv = theLattice->GetIVEnergy(i);		// IV phonon energy

    // Is electron energy above IV Threshold?
    if (eTrk <= Emin_iv) {
        IVprob.push_back(0.); 
        continue;		
    }
      
    G4double scale = 0.;		// IV rate constants
    G4double Efunc = 0.;		// Energy dependence of rates
    G4double ivrate = 0.;		// IV rate
      
    G4double D_iv = theLattice->GetIVDeform(i);		// IV deformation potential
    G4double nVal = theLattice->GetIVNValleys(i);		// # final valleys
    G4double ivorder = theLattice->GetIVOrder(i);		// IV rate order
      
    // 0th order IV rate
    if (ivorder==0) {  
        scale = m_DOS3half * nVal * D_iv*D_iv/ 
            (sqrt(2)*pi*hbar_sq*density*Emin_iv);
        Efunc = energyFunc0th(eTrk-Emin_iv);
    }
      
    // 1st order IV rate
    if (ivorder==1) {  
        G4double qmax = kmag*(1+sqrt(1-Emin_iv/eTrk));		// max phonon momentum
        G4double qmin = kmag*(1-sqrt(1-Emin_iv/eTrk));		// min phonon momentum
        scale = m_DOS3half*m_DOS * nVal * D_iv*D_iv
        /(2*pi*hbar_Planck*density*m_electron*sqrt(m_electron)*Emin_iv*kmag);
        Efunc = energyFunc1st(qmax,qmin);
    }  

    ivrate = scale * Efunc;
    IVprob.push_back(ivrate);
    totalIVRate += ivrate;
  }

  if (verboseLevel>2) G4cout << "IV phonons  " 
      << totalIVRate/hertz << " Hz" << G4endl;
    
  return totalIVRate;
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

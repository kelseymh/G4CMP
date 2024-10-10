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
  //nValley = theLattice->NumberOfValleys()-1;		// From symmetry

  m_DOS = theLattice->GetElectronDOSMass();
  m_DOS3half = sqrt(m_DOS*m_DOS*m_DOS);


}


// Scattering rate is computed from matrix elements

G4double G4CMPInterValleyRate::Rate(const G4Track& aTrack) const {
  const_cast<G4CMPInterValleyRate*>(this)->LoadDataForTrack(&aTrack);
  // Initialize numerical buffers
  eTrk = GetKineticEnergy(aTrack);
  //eTrk = (-1 + sqrt(1 + 4*alpha*GetKineticEnergy(aTrack)))/2/alpha;
    
  ivalley = GetValleyIndex(aTrack);
  ptrk = GetLocalMomentum(aTrack);
  ktrk = theLattice->MapPtoK(ivalley, ptrk);
  kHV = theLattice->EllipsoidalToSphericalTranformation(ivalley, ktrk);
  kmag = kHV.mag();


    
  if (verboseLevel>1)
    G4cout << "G4CMPInterValleyRate eTrk " << eTrk/eV << " eV" << G4endl;

  G4double orate = opticalRate();
  if (verboseLevel>2) G4cout << "IV phonons  " << orate/hertz << " Hz" << G4endl;
 
  G4double nrate = scatterRate();
  if (verboseLevel>2) G4cout << "IV neutrals " << nrate/hertz << " Hz" << G4endl;

  G4double rate =  orate;
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

  IVprob = {};
  std::vector<G4int> nvalleys {7,7};
  std::vector<G4int> orders {0,0};
  std::vector<G4double> ivdeforms {3e8*eV/cm,0.2e8*eV/cm};
    
//   std::vector<G4int> nvalleys {1,1,1,4,4,4};
//   std::vector<G4int> orders {1,1,0,1,0,0};
//   std::vector<G4double> ivdeforms {4*eV,4*eV,11e8*eV/cm,4*eV,2e8*eV/cm,2e8*eV/cm};
    
  G4double total = 0.;
  G4int N_op = theLattice->GetNIVDeform();
  for (G4int i = 0; i<N_op; i++) {
    G4double Emin_op = theLattice->GetIVEnergy(i);
    if (eTrk <= Emin_op) {
        IVprob.push_back(0.); 
        continue;		// Apply threshold behaviour
    }
      
    G4double oscale = 0.;
    G4double scale = 0.;
    G4double Efunc = 0.;
    G4double orate = 0.;    
      
    G4double D_op = theLattice->GetIVDeform(i);
    G4double nVal = theLattice->GetIVNValleys(i);
    G4double ivorder = theLattice->GetIVOrder(i);
  
    if (ivorder==0) {  
    // no kT dependance at ~ mK temperature for IV rate
    scale = nVal*m_DOS3half / (sqrt(2)*pi*hbar_sq*density);
    oscale = scale * D_op*D_op / Emin_op;
    Efunc = energyFunc(eTrk-Emin_op);	// Energy above threshold
    orate = oscale * Efunc;
    }
      
      
    if (ivorder==1) {  
    G4double qmax = kmag*(1+sqrt(1-Emin_op/eTrk));
    G4double qmin = kmag*(1-sqrt(1-Emin_op/eTrk));
    scale = nVal*m_DOS3half*m_DOS / (2*pi*hbar_Planck*density*m_electron*sqrt(m_electron));
    oscale = scale * D_op*D_op / Emin_op/kmag;
    Efunc = qmax*qmax*qmax*qmax-qmin*qmin*qmin*qmin;	// Energy above threshold
    orate = oscale * Efunc;
    }
      
      G4cout << " oscale[" << i << "] : " << oscale << " Efunc : " << Efunc << " Etrk : "  << eTrk/eV << " phonon rate [" << i << "] : " << orate/hertz << " Hz" << G4endl;
      
    if (verboseLevel>2) {
      G4cout << " oscale[" << i << "] " << oscale << " Efunc " << Efunc
	     << "\n phonon rate [" << i << "] " << orate/hertz << " Hz"
	     << G4endl;
    }
         
    IVprob.push_back(orate);
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

/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPLukeEmissionRate.cc
/// \brief Compute emission rate for Luke-Neganov phonons.
//
// $Id$
//
// 20170815  Drop call to LoadDataForTrack(); now handled in process.
// 20170913  Check for electric field; compute "rate" to get up to Vsound
// 20170917  Add interface for threshold identification
// 20240507  Use proper mass and sound of speed for rate.

#include "G4CMPLukeEmissionRate.hh"
#include "G4CMPGeometryUtils.hh"
#include "G4CMPUtils.hh"
#include "G4LatticePhysical.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include <math.h>

#include "G4LatticeLogical.hh"


// Scattering rate is computed from electric field

G4double G4CMPLukeEmissionRate::Rate(const G4Track& aTrack) const {
  // Sanity check -- IsApplicable() should protect against this
  if (!G4CMP::IsChargeCarrier(aTrack)) {
    G4Exception("G4CMPLukeEmissionRate::Rate", "Luke001", EventMustBeAborted, 
		("Invalid particle "+aTrack.GetDefinition()->GetParticleName()).c_str());
    return 0.;
  }

  G4double kmag = 0.; G4double l0 = 0.; G4double mass = 0.; G4double uSound = 0.; G4double massDOS= 0.; G4double alpha = 0; G4double mq = 0; G4double ml = 0; G4double mt = 0; G4double kl = 0; G4double kt = 0; G4ThreeVector kvec; 
  if (G4CMP::IsElectron(aTrack)) {
    kmag = theLattice->MapV_elToK_HV(GetValleyIndex(aTrack),
				     GetLocalVelocityVector(aTrack)).mag();
    kvec = theLattice->MapV_elToK_HV(GetValleyIndex(aTrack),
				     GetLocalVelocityVector(aTrack));
    l0 = theLattice->GetElectronScatter(); 
    massDOS = theLattice->GetElectronDOSMass(); 
    mass = sqrt(electron_mass_c2/c_squared*massDOS); //mass in kSound
    uSound = (2.*theLattice->GetTransverseSoundSpeed() + theLattice->GetSoundSpeed()) / 3.;
    alpha=theLattice->GetAlpha();
  } else if (G4CMP::IsHole(aTrack)) {
    kmag = GetLocalWaveVector(aTrack).mag();
    l0 = theLattice->GetHoleScatter();
    mass = theLattice->GetHoleMass();
    uSound= theLattice->GetSoundSpeed();
  }

  if (verboseLevel > 1) 
    G4cout << "LukeEmissionRate kmag = " << kmag*m << " /m" << G4endl;    
      

  
  G4RotationMatrix massmatrice = theLattice->GetMassTensor();
  ml = massmatrice.xx();
  mt = massmatrice.yy();
  kl = abs(kvec.x());
  kt = sqrt(kvec.y()*kvec.y() + kvec.z()*kvec.z());
  mq=(ml*(kl*kl/kmag/kmag*1/2+kt*kt/kmag/kmag*1/4)+mt*(1-(kl*kl/kmag/kmag*1/2+kt*kt/kmag/kmag*1/4)));
    
  G4double kSound = uSound * sqrt(electron_mass_c2/c_squared*mq) / hbar_Planck;
  //G4double kSound = uSound * mass / hbar_Planck;
    
    
    
        // parabolic, mq case ------------------------------
    
    G4double qmax = 2/(1-2*alpha*mq*uSound*uSound)*(kmag-kSound*sqrt(1+4*alpha*hbar_Planck*hbar_Planck*kmag*kmag/2/electron_mass_c2*c_squared));
    G4double rate = mq*sqrt(mq)/sqrt(massDOS)/massDOS*1/l0/kSound/kSound*uSound/8/kmag*(qmax*qmax*qmax/3*sqrt(1+4*alpha*hbar_Planck*hbar_Planck*kmag*kmag/2/electron_mass_c2*c_squared) - 2*alpha*hbar_Planck*uSound/4*qmax*qmax*qmax*qmax*sqrt(mq/electron_mass_c2*c_squared));
    

    // G4cout << "mq : " << mq/electron_mass_c2*c_squared << G4endl << kvec*m << G4endl << kl*m << G4endl << kt*m  << G4endl; 
    
    
    
    // Non-parabolic, md case ------------------------------
//   G4double qmax = 2/(1-2*alpha*massDOS*uSound*uSound)*(kmag-kSound*sqrt(1+4*alpha*hbar_Planck*hbar_Planck*kmag*kmag/2/electron_mass_c2*c_squared));
//   G4double rate = 1/l0/kSound/kSound*uSound/8/kmag*(qmax*qmax*qmax/3*sqrt(1+4*alpha*hbar_Planck*hbar_Planck*kmag*kmag/2/electron_mass_c2*c_squared) - 2*alpha*hbar_Planck*uSound/4*qmax*qmax*qmax*qmax*sqrt(massDOS/electron_mass_c2*c_squared));
    
    
    
    
    // Print tests ------------------------------
  /* G4cout << "LukeEmissionRate kmag = " << kmag*m << " /m" << G4endl << "rate : " << rate*s << G4endl << "qmax : " << qmax*m << G4endl  <<  "kseuil : " << kSound*sqrt(1+4*alpha*hbar_Planck*hbar_Planck*kmag*kmag/2/electron_mass_c2*c_squared)*m  << G4endl; */  


  // Time step corresponding to Mach number (avg. time between radiations)
  //return (kmag > kSound) ? 1./ChargeCarrierTimeStep(kmag/kSound, l0, uSound) : 0.;
    return (kmag > kSound*sqrt(1+4*alpha*hbar_Planck*hbar_Planck*kmag*kmag/2/electron_mass_c2*c_squared)) ? rate : 0.;
}


// Energy threshold occurs at sound speed in material

G4double G4CMPLukeEmissionRate::Threshold(G4double Eabove) const {
  const G4Track* trk = GetCurrentTrack();	// For convenience below

  G4ThreeVector vtrk = GetLocalVelocityVector(trk);
  G4double vsound = theLattice->GetSoundSpeed();
  G4ThreeVector v_el = vsound * vtrk.unit();
  G4double Esound = theLattice->MapV_elToEkin(GetValleyIndex(trk), v_el);

  if (verboseLevel>1) {
    G4cout << "G4CMPLukeEmissionRate::Threshold vtrk " << vtrk.mag()/(m/s)
	   << " vsound " << vsound/(m/s) << " m/s Esound " << Esound/eV
	   << " eV" << G4endl;
  }

  // Thresholds or pseudothresholds at multiples of Esound
  const G4double eStep = 25.;
  G4double ratio = Eabove/Esound;
  if (ratio > 1.) {
    ratio += (std::fmod(ratio,eStep)<0.95) ? 0. : eStep;  // Avoid Zeno paradox
    ratio = std::ceil(ratio/eStep) * eStep;

    if (verboseLevel>2) G4cout << " scaling by " << ratio << G4endl;
    
    Esound *= ratio;
  }

  return (Eabove < Esound) ? Esound : 0.;	// No thresholds above sound
}

//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file library/src/G4eLukeScattering.cc
/// \brief Implementation of the G4eLukeScattering class
//
// $Id$
//

#include "G4CMPeLukeScattering.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"

#include "G4PhononLong.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"

#include "G4FieldManager.hh"
#include "G4Step.hh"
#include "G4VParticleChange.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "G4TransportationManager.hh"
#include "G4CMPValleyTrackMap.hh"
#include "G4Field.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include <fstream>
#include <iostream>


G4CMPeLukeScattering::G4CMPeLukeScattering(G4VProcess* sLim)
 :G4VPhononProcess("eLukeScattering"), stepLimiter(sLim),
  velLong(5324.2077*m/s), l0(257e-6*m), me(electron_mass_c2/c_squared),
  mc(.1185*me), ksound(mc*velLong/hbar_Planck) {
  G4cout << "G4CMPeLukeScattering::Constructor: ksound = "
	 << ksound*m << " /m" << G4endl;

  if (verboseLevel > 1) {
    G4cout << GetProcessName() << " is created." << G4endl;
  }
}

G4CMPeLukeScattering::~G4CMPeLukeScattering() {;}


G4double G4CMPeLukeScattering::GetMeanFreePath(const G4Track& aTrack,
					       G4double,
					       G4ForceCondition* condition) {
  G4RotationMatrix trix;
  int valley = G4CMPValleyTrackMap::GetInstance()->GetValley(aTrack);
    
  switch (valley) {
  case 1:
    trix = G4RotationMatrix(
	     G4Rep3x3( 1.0/sqrt(3.0),  1.0/sqrt(3.0),  1.0/sqrt(3.0),
		      -1.0/sqrt(2.0),  1.0/sqrt(2.0),  0.0,
		      -1.0/sqrt(6.0), -1.0/sqrt(6.0),  sqrt(2.0/3.0) ));
    break;
  case 2:
    trix = G4RotationMatrix(
	     G4Rep3x3(-1.0/sqrt(3.0),  1.0/sqrt(3.0),  1.0/sqrt(3.0),
		      -1.0/sqrt(2.0), -1.0/sqrt(2.0),  0.0,
		       1.0/sqrt(6.0), -1.0/sqrt(6.0),  sqrt(2.0/3.0) ));
    break;
  case 3:
    trix = G4RotationMatrix(
	     G4Rep3x3(-1.0/sqrt(3.0), -1.0/sqrt(3.0),  1.0/sqrt(3.0),
		       1.0/sqrt(2.0), -1.0/sqrt(2.0),  0.0,
		       1.0/sqrt(6.0),  1.0/sqrt(6.0),  sqrt(2.0/3.0) ));
    break;
  case 4:
    trix = G4RotationMatrix(
	     G4Rep3x3( 1.0/sqrt(3.0), -1.0/sqrt(3.0),  1.0/sqrt(3.0),
		       1.0/sqrt(2.0),  1.0/sqrt(2.0),  0.0,
		      -1.0/sqrt(6.0),  1.0/sqrt(6.0),  sqrt(2.0/3.0) ));
    break;
  }

  valleyToNormal= G4AffineTransform(trix);
  normalToValley= G4AffineTransform(trix).Inverse();    
  T = G4ThreeVector(sqrt(.118/1.588), sqrt(.118/.081), sqrt(.118/.081));


  G4ThreeVector k = aTrack.GetMomentum()/hbarc;
  G4ThreeVector v = mInv*k*hbar_Planck;
  G4ThreeVector k_valley = normalToValley.TransformPoint(k);
  G4ThreeVector k_HV =  G4ThreeVector(k_valley[0]*T[0],
				      k_valley[1]*T[1],
				      k_valley[2]*T[2]);
  
  *condition = Forced;
  
  G4double kmag = k_HV.mag();
  
  if (kmag<=ksound) return DBL_MAX;
  
  G4double tau0 =  1.0 / (
			  velLong / (3*l0)
			  * (kmag / ksound) * (kmag/ksound)
			  * ((1- ksound /kmag))
			  * ((1- ksound /kmag))
			  * ((1- ksound /kmag))
			  );
  
  G4double mfp0 = tau0*v.mag();
  return mfp0;
}

G4VParticleChange* G4CMPeLukeScattering::PostStepDoIt(const G4Track& aTrack,
						      const G4Step& aStep) {
  aParticleChange.Initialize(aTrack); 
  G4StepPoint* postStepPoint = aTrack.GetStep()->GetPostStepPoint();
  
  G4ThreeVector k = postStepPoint->GetMomentum()/hbarc;
  G4ThreeVector k_valley = normalToValley.TransformPoint(k);
  G4ThreeVector k_HV =  G4ThreeVector(k_valley[0]*T[0],
				      k_valley[1]*T[1],
				      k_valley[2]*T[2]);
  
  G4double kmag = k_HV.mag();
  
  //Do nothing other than re-calculate mfp when step limit reached or leaving
  //volume
  if ( (postStepPoint->GetProcessDefinedStep()==stepLimiter)
       || (postStepPoint->GetStepStatus()==fGeomBoundary)
       || (kmag <= ksound) ) {
    return &aParticleChange;
  }

  G4double theta_phonon = MakeTheta(kmag, ksound);
  G4double theta_charge = 
    acos( (kmag*kmag - 2*ksound*(kmag*cos(theta_phonon)
				 - ksound)
	   - 2 * (kmag*cos(theta_phonon) - ksound)
	   * (kmag*cos(theta_phonon) - ksound) )/
	  kmag/ (sqrt(kmag*kmag - 4*ksound
		      * (kmag*cos(theta_phonon) - ksound) ) ) );
  
  G4double q = 2*(kmag*cos(theta_phonon)-ksound);
  //k_HV.setMag(sqrt(k_HV.mag2() + q*q - 2*kmag*q*cos(theta_phonon)));
  k_HV.setMag(sqrt(kmag*kmag-2*mc*q*velLong/hbar_Planck));
  G4ThreeVector kdir = k_HV.unit();
  k_HV.rotate(kdir.orthogonal(), theta_charge);
  
  G4double phi_charge =  G4UniformRand()*2*pi;
  k_HV.rotate(kdir, phi_charge);
  
  G4ThreeVector p_new = hbar_Planck*k_HV;
  p_new[0] /= T[0];
  p_new[1] /= T[1];
  p_new[2] /= T[2];
  valleyToNormal.ApplyPointTransform(p_new);
  
  G4ThreeVector phononq = q*k.unit().rotate(k.unit(), theta_phonon);
  //   std::ofstream charge("theta_charge", std::ios::app);
  std::ofstream phonon("theta_phonon", std::ios::app);
  //   charge << theta_charge << G4endl;
  phonon << phononq.getZ()/phononq.mag()<< G4endl;
  //   charge.close();
  phonon.close();

  std::ofstream espec("energy_spectrum_phonon", std::ios::app);
  espec << (k_HV.mag2()*hbar_Planck/2/mc)/eV*1000<< G4endl;
  espec.close();
  
  aParticleChange.ProposeMomentumDirection(p_new.unit());
  aParticleChange.ProposeEnergy(p_new.mag2()/2/mc);
  ResetNumberOfInteractionLengthLeft();
  return &aParticleChange;
}


G4double G4CMPeLukeScattering::MakeTheta(G4double& k, G4double& ks) {
  //Analytical method to compute theta
  G4double u = G4UniformRand();
  G4double base = ( -(u-1)+3*(u-1)*(ks/k)
		    -3*(u-1)*(ks/k)*(ks/k)
		    +(u-1)*(ks/k)*(ks/k)*(ks/k) );
  G4double exponent = 1.0/3.0;
  G4double operand = ks/k+pow(base, exponent);
  
  if(operand>1.0) {
        // G4cout<<"\nTruncating operand from"<<operand<<" to 1.0";
        operand=1.0;
    }
  if(operand<0.0) G4cout<<"\noperand out of range: operand = "<<operand;
  
  G4double theta = acos(operand);
  
  if(acos(ks/k)<theta) G4cout<<"\n  THETA OUT OF BOUNDS!!! (theta>acos(ks/k)) theta:"<<theta<<" acos(ks/k):"<<acos(ks/k);
  if(pi/2<=theta) G4cout<<"\n THETA OUT OF BOUNDS!!! (pi/2 < theta)";
  
  return theta;
}

G4double G4CMPeLukeScattering::MakePhi(G4double& k,G4double& ks, G4double& theta) {
  G4double phi=acos((k*k - 2*ks*(k*cos(theta)-ks)-2*(k*cos(theta)-ks)*(k*cos(theta)-ks))/(k*sqrt(k*k-4*ks*(k*cos(theta)-ks))));
  return phi;
}

G4bool G4CMPeLukeScattering::IsApplicable(const G4ParticleDefinition& aPD) {
  return((&aPD==G4CMPDriftElectron::Definition()));
}

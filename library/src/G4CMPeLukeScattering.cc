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
 :G4VPhononProcess("eLukeScattering")
{
  velLong=5324.2077*m/s;
  l0Hole = 257e-6*m;
  massFreeElectron=9.1e-31*kg;
  massHole=.12*massFreeElectron;

  //ksound_Hole=5512124.0235/m; //velLong*massHole / hbar;//1.6306e7/m;
  ksound_Hole = velLong*massHole/hbar_Planck;
  G4cout<<"\nG4CMPeLukeScattering::Constructor:ksound_Hole = "
	<<ksound_Hole*m<<" /m"<<G4endl;

	stepLimiter = sLim;
  /*
    normalToValley= G4AffineTransform
    (G4RotationMatrix(G4ThreeVector(1.0,0,0), 0.95532));
  valleyToNormal= G4AffineTransform
    (G4RotationMatrix(G4ThreeVector(1.0,0,0), 0.95532)).Inverse();
  */

  /*
  HVMatrix = G4RotationMatrix(G4ThreeVector(1.2172, 0     ,     0),
				      G4ThreeVector(0     , 1.2172,     0),
				      G4ThreeVector(0     , 0     ,0.2756)
				      );
  HVTransform = G4AffineTransform(HVMatrix);
  */

  if(verboseLevel>1){
    G4cout<<GetProcessName()<<" is created "<<G4endl;
  }
}

G4CMPeLukeScattering::~G4CMPeLukeScattering()
{ ; }

G4double G4CMPeLukeScattering::GetMeanFreePath(const G4Track& aTrack, G4double, G4ForceCondition* condition)
{
  //  G4cout << "GetMeanFreePath" << G4endl;
  G4double me=electron_mass_c2/c_squared;
  G4double mc = .118*me;

  G4RotationMatrix trix;
  int valley = G4CMPValleyTrackMap::GetInstance()->GetValley(aTrack);
    
  switch(valley){
      case 1:
	  trix = G4RotationMatrix(G4Rep3x3( 1.0/sqrt(3.0),  1.0/sqrt(3.0),  1.0/sqrt(3.0), 
			  -1.0/sqrt(2.0),  1.0/sqrt(2.0),  0.0, 
			  -1.0/sqrt(6.0), -1.0/sqrt(6.0),  sqrt(2.0/3.0) ));
          break;
      case 2:
	  trix = G4RotationMatrix(G4Rep3x3(-1.0/sqrt(3.0),  1.0/sqrt(3.0),  1.0/sqrt(3.0), 
		 -1.0/sqrt(2.0), -1.0/sqrt(2.0),  0.0, 
		  1.0/sqrt(6.0), -1.0/sqrt(6.0),  sqrt(2.0/3.0) ));
          break;
      case 3:
	  trix = G4RotationMatrix(G4Rep3x3(-1.0/sqrt(3.0), -1.0/sqrt(3.0),  1.0/sqrt(3.0), 
		  1.0/sqrt(2.0), -1.0/sqrt(2.0),  0.0, 
		  1.0/sqrt(6.0),  1.0/sqrt(6.0),  sqrt(2.0/3.0) ));
          break;
      case 4:
	  trix = G4RotationMatrix(G4Rep3x3( 1.0/sqrt(3.0), -1.0/sqrt(3.0),  1.0/sqrt(3.0), 
		  1.0/sqrt(2.0),  1.0/sqrt(2.0),  0.0, 
		 -1.0/sqrt(6.0),  1.0/sqrt(6.0),  sqrt(2.0/3.0) ));
          break;
  }
  
  valleyToNormal= G4AffineTransform(trix);
  normalToValley= G4AffineTransform(trix).Inverse();    
  G4ThreeVector T = G4ThreeVector(sqrt(.118/1.588), sqrt(.118/.081), sqrt(.118/.081));
  G4FieldManager* fMan = 
	G4TransportationManager::GetTransportationManager()->GetFieldManager();
    const G4Field* field = fMan->GetDetectorField();

    G4ThreeVector posVec = aTrack.GetPosition();
    G4double position[4] = {posVec[0],posVec[1],posVec[2],0};
    G4double fieldVal[6];

    field->GetFieldValue(position,fieldVal);
    G4ThreeVector Efield = G4ThreeVector(fieldVal[3], fieldVal[4], fieldVal[5]);
    G4ThreeVector Efield_valley = normalToValley.TransformPoint(Efield);
    G4ThreeVector Efield_HV = G4ThreeVector( Efield_valley[0]*T[0], 
					     Efield_valley[1]*T[1], 
					     Efield_valley[2]*T[2]);
    
    G4RotationMatrix mInv = 
		trix.inverse()*G4Rep3x3(1/1.588/me,   0.0    , 0.0,
							0.0     , 1/.081/me, 0.0, 
							0.0     ,   0.0    , 1/.081/me)
							*trix;

  G4ThreeVector k = aTrack.GetMomentum()/hbarc;
  G4ThreeVector v = mInv*k*hbar_Planck;
  G4ThreeVector k_valley = normalToValley.TransformPoint(k);
  G4ThreeVector k_HV =  G4ThreeVector(k_valley[0]*T[0], 
					k_valley[1]*T[1], 
					k_valley[2]*T[2]);
					
  *condition = Forced;

  G4double kmag = k_HV.mag();
					
  if(kmag<=ksound_Hole) 
  {
      return DBL_MAX;
  }

  G4double tau0 =  1.0 / (
		       velLong / (3*l0Hole)
		       * (kmag / ksound_Hole) * (kmag/ksound_Hole)
		       * ((1- ksound_Hole /kmag))
		       * ((1- ksound_Hole /kmag))
		       * ((1- ksound_Hole /kmag))
		       );
	
  G4double mfp0 = tau0*v.mag();
  return mfp0;
  G4cout << mfp0/m << G4endl;
  
  G4ThreeVector k1, k1_HV, x1, v1;
  G4ThreeVector v0 = v;
  G4ThreeVector x0 = posVec;
  G4double mfp1, tau1;
  G4double dmfp = 0;
  for (int i=0; i<100; i++)
  {
  x1 = x0 + mfp0*v0/v0.mag();
  k1_HV =  k_HV + electron_charge/hbar_Planck*Efield_HV*mfp0/v0.mag()/100;
  v1 = v0 + mInv*Efield*electron_charge*mfp0/v0.mag()/100;
  
  //G4cout << (k1_HV.mag() - k_HV.mag())*m << G4endl;
  
  tau1 =  1.0 / (
		       velLong / (3*l0Hole)
		       * (k1_HV.mag()/ ksound_Hole)
		       * ((1- ksound_Hole /k1_HV.mag()))
		       * ((1- ksound_Hole /k1_HV.mag()))
		       * ((1- ksound_Hole /k1_HV.mag()))
		       );

    //G4cout << "mfp0 = " << mfp0/m << G4endl;
    //G4cout << "mfp1 = " << tau1*v1.mag()/m << G4endl;
    //G4cout << "mfp2 = " << tau1*v1.mag()/m  + mfp0/10/m<< G4endl;
    dmfp += mfp0/100;
    mfp0 = mfp0/100 + tau1*v1.mag();
    x0 = x1;
    k_HV = k1_HV;
    v0 = v1;
  }
  //dmfp += tau1*v1.mag();
//   for (int i=0; i<100; i++)
//   {
//   //x1 = x0 + mfp0*v0/v0.mag();
//   k1_HV =  k_HV + electron_charge/hbar_Planck*Efield_HV*mfp0/v0.mag()/100;
//   v1 = v0 + mInv*Efield*electron_charge*mfp0/v0.mag()/100;
//   
//   G4cout << (k1_HV.mag() - k_HV.mag())*m << G4endl;
//   
//   tau1 =  1.0 / (
// 		       velLong / (3*l0Hole)
// 		       * (k1_HV.mag()/ ksound_Hole)
// 		       * ((1- ksound_Hole /k1_HV.mag()))
// 		       * ((1- ksound_Hole /k1_HV.mag()))
// 		       * ((1- ksound_Hole /k1_HV.mag()))
// 		       );
// 
//     G4cout << "mfp0 = " << mfp0/m << G4endl;
//     G4cout << "mfp1 = " << tau1*v1.mag()/m << G4endl;
//     G4cout << "mfp2 = " << tau1*v1.mag()/m  + mfp0/10/m<< G4endl;
//     mfp0 = mfp0/100 + tau1*v1.mag();
//   }
 G4cout << mfp0/m << G4endl; 
  return dmfp;
  
  
}

G4VParticleChange* G4CMPeLukeScattering::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
  G4double me=electron_mass_c2/c_squared;
  G4double mc = .118*me;
  aParticleChange.Initialize(aTrack); 
  G4RotationMatrix trix;
  int valley = G4CMPValleyTrackMap::GetInstance()->GetValley(aTrack);
    
  switch(valley){
      case 1:
	  trix = G4RotationMatrix(G4Rep3x3( 1.0/sqrt(3.0),  1.0/sqrt(3.0),  1.0/sqrt(3.0), 
			  -1.0/sqrt(2.0),  1.0/sqrt(2.0),  0.0, 
			  -1.0/sqrt(6.0), -1.0/sqrt(6.0),  sqrt(2.0/3.0) ));
          break;
      case 2:
	  trix = G4RotationMatrix(G4Rep3x3(-1.0/sqrt(3.0),  1.0/sqrt(3.0),  1.0/sqrt(3.0), 
		 -1.0/sqrt(2.0), -1.0/sqrt(2.0),  0.0, 
		  1.0/sqrt(6.0), -1.0/sqrt(6.0),  sqrt(2.0/3.0) ));
          break;
      case 3:
	  trix = G4RotationMatrix(G4Rep3x3(-1.0/sqrt(3.0), -1.0/sqrt(3.0),  1.0/sqrt(3.0), 
		  1.0/sqrt(2.0), -1.0/sqrt(2.0),  0.0, 
		  1.0/sqrt(6.0),  1.0/sqrt(6.0),  sqrt(2.0/3.0) ));
          break;
      case 4:
	  trix = G4RotationMatrix(G4Rep3x3( 1.0/sqrt(3.0), -1.0/sqrt(3.0),  1.0/sqrt(3.0), 
		  1.0/sqrt(2.0),  1.0/sqrt(2.0),  0.0, 
		 -1.0/sqrt(6.0),  1.0/sqrt(6.0),  sqrt(2.0/3.0) ));
          break;
  }
  
  
  G4StepPoint* postStepPoint = aTrack.GetStep()->GetPostStepPoint();
  valleyToNormal= G4AffineTransform(trix);
  normalToValley= G4AffineTransform(trix).Inverse();
    
    G4RotationMatrix mInv = 
		trix.inverse()*G4Rep3x3(1/1.588/me,   0.0    , 0.0,
							0.0     , 1/.081/me, 0.0, 
							0.0     ,   0.0    , 1/.081/me)
							*trix;
  
  G4ThreeVector T = G4ThreeVector(sqrt(.12/1.588), sqrt(.12/.081), sqrt(.12/.081));
  G4ThreeVector k = postStepPoint->GetMomentum()/hbarc;
  G4ThreeVector k_valley = normalToValley.TransformPoint(k);
  G4ThreeVector k_HV =  G4ThreeVector(k_valley[0]*T[0], 
					k_valley[1]*T[1], 
					k_valley[2]*T[2]);
					
  G4double kmag = k_HV.mag(); 
  
  //Do nothing other than re-calculate mfp when step limit reached or leaving 
  //volume
  if((postStepPoint->GetProcessDefinedStep()==stepLimiter)
||(postStepPoint->GetStepStatus()==fGeomBoundary)
       || (kmag <= ksound_Hole)
     )
    {     
	return &aParticleChange;
    }  

  //G4cout << "Luke Emission1" << G4endl;
  G4double theta_phonon= MakeTheta(kmag, ksound_Hole);

  G4double theta_charge=acos( 
			     (kmag*kmag - 2*ksound_Hole
			      *(kmag*cos(theta_phonon) - ksound_Hole) 
			      - 2 * (kmag*cos(theta_phonon) - ksound_Hole)
			      * (kmag*cos(theta_phonon) - ksound_Hole)
			      ) 
			     / kmag
			     / (sqrt(
				     kmag*kmag - 4*ksound_Hole
				     *(kmag*cos(theta_phonon) - ksound_Hole)
				     )) 
			     );

  G4double q = 2*(kmag*cos(theta_phonon)-ksound_Hole);
  k_HV.setMag(sqrt(k_HV.mag2() + q*q - 2*kmag*q*cos(theta_phonon)));
  //k_HV.setMag(sqrt(kmag*kmag-2*mc*q*velLong/hbar_Planck));

  G4ThreeVector kdir = k_HV.unit();
  //G4cout << "kmag = " << k_HV.mag() << G4endl;
  k_HV.rotate(kdir.orthogonal(), theta_charge);
  //G4cout << "kmag = " << k_HV.mag() << G4endl;

  G4double phi_charge =  G4UniformRand()*2*pi;
  k_HV.rotate(kdir, phi_charge);
  //G4cout << "kmag = " << k_HV.mag() << G4endl;
  
G4ThreeVector p_new = hbar_Planck*k_HV;
p_new[0] /= T[0];
p_new[1] /= T[1];
p_new[2] /= T[2];
valleyToNormal.ApplyPointTransform(p_new);

/*
  std::ofstream epositions;
  epositions.open("qenergies.txt", std::ofstream::app);
files
  epositions << aTrack.GetKineticEnergy()/eV - Energy << "\n";
  epositions.close();
  */


//    std::ofstream energy;
//    energy.open("energy", std::ofstream::app);
//  
//    energy<< aTrack.GetKineticEnergy()/eV << "\n";
//    energy.close();
   
    std::ofstream qwave;
    qwave.open("q", std::ofstream::app);
    qwave << hbar_Planck*velLong*q/eV << "\n";
    qwave.close();
//   
//   std::ofstream kwave;
//   kwave.open("k", std::ofstream::app);
// 
//   kwave << oldE/eV - Energy/eV << "\n";
//   kwave.close();
/*
G4cout << "Energy before   : " << aTrack.GetKineticEnergy()/eV << G4endl;
G4cout << "Energy before   : " << k.mag2()*hbar_Planck*hbar_Planck/2/mc/eV<< G4endl;
G4cout << "Energy after    : " << p_new.mag2()/2/mc/eV << G4endl;
G4cout << "Energy of Phonon: " << velLong*hbar_Planck*q/eV << G4endl;*/
  
  
  //G4cout << postStepPoint->GetVelocity()*s/m << G4endl;
  //G4cout << v_new.mag()*s/m << G4endl;
  aParticleChange.ProposeMomentumDirection(p_new.unit());
  //aParticleChange.ProposeVelocity(v_new.mag()); 
  //aParticleChange.ProposeEnergy(Energy);
  //aParticleChange.ProposeEnergy(p_new*(mInv*p_new)/2);
  aParticleChange.ProposeEnergy(p_new.mag2()/2/mc);
  //aParticleChange.ProposeEnergy(aTrack.GetKineticEnergy()-velLong*hbar_Planck*q);
  //aParticleChange.ProposeEnergy(T -velLong*hbar_Planck*q);
  ResetNumberOfInteractionLengthLeft();
  return &aParticleChange;
}

G4double G4CMPeLukeScattering::MakeTheta(G4double& k, G4double& ks){
   /*  G4double u = G4UniformRand();
 
  //Analytical method to compute theta
  

  double base = -(u-1)+3*(u-1)*(ks/k)-3*(u-1)*(ks/k)*(ks/k)+(u-1)*(ks/k)*(ks/k)*(ks/k);
  double exponent =1.0/3.0;
  G4double operand = ks/k+pow(base, exponent);   

  if(operand>1.0) {
    // G4cout<<"\nTruncating operand from"<<operand<<" to 1.0";
    operand=1.0;
  }
  if(operand<0.0) G4cout<<"\noperand out of range: operand = "<<operand;
  G4double theta=acos(operand);
  //  G4cout<<"\n"<<theta;

  if(acos(ks/k)<theta) G4cout<<"\n  THETA OUT OF BOUNDS!!! (theta>acos(ks/k)) theta:"<<theta<<" acos(ks/k):"<<acos(ks/k);
  if(PI/2<=theta) G4cout<<"\n THETA OUT OF BOUNDS!!! (pi/2 < theta)";
      
      */
  //Rejection method for determining theta
  G4double theta;
  G4double pValue=2;
  G4double pDensity = 1;

  G4double thetaMax;
  G4double pValMax;
  thetaMax=acos(ks/k);
  pValMax=(1-(ks/k))*(1-(ks/k));
  bool first = false;

  while(pValue>pDensity){
    pValue=G4UniformRand();//*pValMax;//  *(1+2*ks/k+(ks/k)*(ks/k));   
    theta=G4UniformRand()*thetaMax;
    pDensity = (cos(theta)-(ks/k))*(cos(theta)-(ks/k))*sin(theta);//need to multiply by unit 's' to make it dimensionless
    if(pDensity>pValMax) G4cout<<"\nLukeScattering::PostStepDoIt: Error: pDensity should never exceed pValMax "<<pValMax;
    
    
    //    if(!first){
    //G4cout<<"\nG4CMPLukeScattering::MakeTheta: pDensity calculated as: 
//"<<pDensity;
  //   G4cout<<"\n\tpValue: "<<pValue;
    // G4cout<<"\n\ttheta: "<<theta/rad;
    // G4cout<<"\n\tk: "<<k*m;
   //  G4cout<<"\n\tks: "<<ks*m;
   //  G4cout<<endl;
    //}
    //first=false;
    
  }
  return theta;
}

G4double G4CMPeLukeScattering::MakePhi(G4double& k,G4double& ks, G4double& theta){


  G4double phi=acos((k*k - 2*ks*(k*cos(theta)-ks)-2*(k*cos(theta)-ks)*(k*cos(theta)-ks))/(k*sqrt(k*k-4*ks*(k*cos(theta)-ks))));

  return phi;
}

G4bool G4CMPeLukeScattering::IsApplicable(const G4ParticleDefinition& aPD)
{
  
return((&aPD==G4CMPDriftElectron::Definition()));
}

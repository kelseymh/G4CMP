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
/// \file library/src/G4PhononReflection.cc
/// \brief Implementation of the G4PhononReflection class
//
// This process handles the interaction of phonons with
// boundaries. Implementation of this class is highly 
// geometry dependent.Currently, phonons are killed when
// they reach a boundary. If the other side of the 
// boundary was Al, a hit is registered.
//  
// $Id$
//
// 20131115  Throw exception if track's polarization state is invalid.
// 20140103  Move charge version of code here, still commented out
// 20140331  Add required process subtype code

#include "G4PhononReflection.hh"
#include "G4ExceptionSeverity.hh"
#include "G4GeometryTolerance.hh"
#include "G4LatticePhysical.hh"
#include "G4PhononLong.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononTransSlow.hh"
#include "G4PhysicalConstants.hh"
#include "G4Navigator.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4TransportationManager.hh"
#include "G4VParticleChange.hh"
#include "Randomize.hh"


G4PhononReflection::G4PhononReflection(const G4String& aName)
  : G4VPhononProcess(aName, fPhononReflection),
    kCarTolerance(G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()) {;}

G4PhononReflection::~G4PhononReflection() {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Always return DBL_MAX and Forced. This ensures that the process is
// called at the end of every step. In PostStepDoIt the process
// decides whether the step encountered a volume boundary and a
// reflection should be applied

G4double G4PhononReflection::GetMeanFreePath(const G4Track&, G4double,
					     G4ForceCondition* condition) {
  *condition = Forced;
  return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// This process handles the interaction of phonons with
// boundaries. Implementation of this class is highly geometry
// dependent.Currently, phonons are killed when they reach a
// boundary. If the other side of the boundary was Al, a hit is
// registered.
  
G4VParticleChange* G4PhononReflection::PostStepDoIt(const G4Track& aTrack,
						    const G4Step& aStep) { 
  aParticleChange.Initialize(aTrack);
   
  //Check if current step is limited by a volume boundary
  G4StepPoint* postStepPoint = aStep.GetPostStepPoint();
  if (postStepPoint->GetStepStatus()!=fGeomBoundary) {
    //make sure that correct phonon velocity is used after the step
    int pol = GetPolarization(aTrack);
    if (pol < 0 || pol > 2) {
      G4Exception("G4PhononReflection::PostStepDoIt","Phonon001",
		  EventMustBeAborted, "Track is not a phonon");
      return &aParticleChange;		// NOTE: Will never get here
    }

    // FIXME:  This should be using wave-vector, shouldn't it?
    G4double vg = theLattice->MapKtoV(pol, aTrack.GetMomentumDirection());
    
    //Since step was not a volume boundary, just set correct phonon velocity and return
    aParticleChange.ProposeVelocity(vg);
    return &aParticleChange;
  }
  
  // do nothing but return is the step is too short
  // This is to allow actual reflection where after
  // the first boundary crossing a second, infinitesimal
  // step occurs crossing back into the original volume
  if (aTrack.GetStepLength()<=kCarTolerance/2) { 
    return &aParticleChange;
  }

  /*** THIS IS WHERE REFLECTION/TRANSMISSION CALCULATIONS BELONG ***

  const G4DynamicParticle* theDP = aTrack.GetDynamicParticle();
  G4ThreeVector incidentDirection = theDP->GetMomentumDirection();
  G4Navigator* theNavigator = 
    G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();

  G4bool valid = true;
  G4ThreeVector localNormal = theNavigator->GetLocalExitNormal(&valid);
  if (valid) { localNormal = -localNormal; }
  G4ThreeVector globalNormal = 
    theNavigator->GetLocalToGlobalTransform().TransformAxis(localNormal);
  
  //Set specular scattering probability SSP//
  G4double SSP = 0;

  if (G4UniformRand()<SSP){
    if (incidentDirection*globalNormal>0.0) {
      // this should not happen but .......
      globalNormal = - globalNormal;
    }

    G4double PdotN = incidentDirection*globalNormal;

    reflectedDirection = incidentDirection - (2.*PdotN)*globalNormal;
    info->setK(reflectedDirection);  //This is just a temporary bug fix, in order to determine a k-vector. Not physical
    //reflectedDirection=G4LatticeManager::mapKtoVDir(aTrack.GetVolume(),2,reflectedDirection);
  } else {
    /////////If scattering is diffuse:///////////
    G4ThreeVector reflectedK;
    
    G4double PdotN = incidentDirection*globalNormal;
    if (PdotN>0.) {
      // this should not happen but .......
      globalNormal = - globalNormal;
      PdotN *= -1.;
    }
    
    reflectedDirection = G4LambertianRand(globalNormal);
    
    info->setK(reflectedDirection); 
    //if(aTrack.GetDefinition()==G4PhononLong::PhononDefinition()) reflectedDirection=G4LatticeManager::mapKtoVDir(aTrack.GetVolume(),0,reflectedDirection);
    //else if(aTrack.GetDefinition()==G4PhononTransSlow::PhononDefinition()) reflectedDirection=G4LatticeManager::mapKtoVDir(aTrack.GetVolume(),1,reflectedDirection);
    //else if(aTrack.GetDefinition()==G4PhononTransFast::PhononDefinition()) reflectedDirection=G4LatticeManager::mapKtoVDir(aTrack.GetVolume(),2,reflectedDirection);
  }
  
  aParticleChange.ProposeMomentumDirection(reflectedDirection.unit());
  
  //   check if phonon is lost to black body radiation
  if (postStepPoint->GetMaterial() != alminum) {
    if (G4UniformRand()<0.001) { 
      aParticleChange.ProposeTrackStatus(fStopAndKill);
    }
  } else {
    // in case the other side is alminum electrode, deposit energy
    //   QPDownconversion con;
    if ((G4UniformRand()<0.02013) && 
	(aTrack.GetKineticEnergy()>347.43e-6*eV)) {
      //       con.downconvert(aTrack, &aParticleChange, reflectedDirection);
      G4double TwoAlGap = 347.43e-6*eV;
      G4double eKin = aTrack.GetKineticEnergy();
      
      if(eKin<4*TwoAlGap){
	aParticleChange.ProposeTrackStatus(fStopAndKill);   
	aParticleChange.ProposeNonIonizingEnergyDeposit(eKin);
      } else {
	aParticleChange.ProposeNonIonizingEnergyDeposit(4*TwoAlGap);
	aParticleChange.ProposeEnergy(eKin-4*TwoAlGap);
	//if(G4UniformRand()>(5.6/15.0)){
	//aParticleChange.ProposeTrackStatus(fStopAndKill);
	//aParticleChange.ProposeEnergy(0);
	//}
      }
    } else if (aTrack.GetKineticEnergy()<347.4e-6*eV){
      aParticleChange.ProposeTrackStatus(fStopAndKill);  
    }
  } else
  ***/
  {
    //*** FOR NOW, PHONONS ARE SIMPLY KILLED AT THE WALL ***
    G4double eKin = aTrack.GetKineticEnergy();     
    aParticleChange.ProposeNonIonizingEnergyDeposit(eKin);
    aParticleChange.ProposeTrackStatus(fStopAndKill);
  }

  return &aParticleChange; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



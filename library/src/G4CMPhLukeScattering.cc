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
/// \file library/src/G4LukeScattering.cc
/// \brief Implementation of the G4LukeScattering class
//
// $Id$
//
// 20140325  Move time-step calculation to G4CMPProcessUtils
// 20140331  Add required process subtype code

#include "G4CMPhLukeScattering.hh"
#include "G4CMPDriftHole.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalVolume.hh"
#include "G4PhysicalConstants.hh"
#include "G4RandomDirection.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4VParticleChange.hh"
#include "G4VPhysicalVolume.hh"
#include "Randomize.hh"
#include <math.h>


G4CMPhLukeScattering::G4CMPhLukeScattering(G4VProcess* stepper)
  : G4CMPVDriftProcess("hLukeScattering", fLukeScattering),
    stepLimiter(stepper) {;}

G4CMPhLukeScattering::~G4CMPhLukeScattering() {;}


G4double 
G4CMPhLukeScattering::GetMeanFreePath(const G4Track& aTrack, G4double,
				      G4ForceCondition* condition) {
  *condition = Forced;
  
  G4StepPoint* stepPoint = aTrack.GetStep()->GetPostStepPoint();
  G4double velocity = stepPoint->GetVelocity();
  G4double kmag = velocity*mc_h / hbar_Planck;
  
  if (kmag<=ksound_h) return DBL_MAX;
  
  // Time step corresponding to Mach number (avg. time between radiations)
  G4double dtau = ChargeCarrierTimeStep(kmag/ksound_h, l0_h);
  
  G4double mfp = dtau * velocity;

#ifdef G4CMP_DEBUG  
  G4cout << "hLuke MFP = " <<  mfp <<  G4endl;
#endif
  return mfp;
}

G4VParticleChange* G4CMPhLukeScattering::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
  aParticleChange.Initialize(aTrack);  
  
  //Do nothing other than re-calculate mfp when step limit reached or leaving volume
  G4StepPoint* postStepPoint = aStep.GetPostStepPoint();
#ifdef G4CMP_DEBUG
  G4cout << GetProcessName() << "PostStepDoIt: Step limited by process "
	 << postStepPoint->GetProcessDefinedStep()->GetProcessName()
	 << G4endl;
#endif
  
  if((postStepPoint->GetProcessDefinedStep()==stepLimiter)
     || (postStepPoint->GetStepStatus()==fGeomBoundary)) {
    return &aParticleChange;
  }
  
  G4double velocity = postStepPoint->GetVelocity();
  //G4double p = postStepPoint->GetMomentum().mag()/c_light;
  //G4cout << "momentum: " <<  p <<  " " << velocity*mc_h << G4endl;
  G4double kmag = velocity*mc_h / hbar_Planck;

  G4double theta_phonon = MakePhononTheta(kmag, ksound_h);
  G4double theta_charge = MakeRecoilTheta(kmag, ksound_h, theta_phonon);
  
  G4double q = 2*(kmag*cos(theta_phonon)-ksound_h);
  G4double T = hbar_Planck*hbar_Planck*kmag*kmag/2/mc_h;
  G4double qEnergy = velLong*hbar_Planck*q;

  G4ThreeVector momentum = G4ThreeVector(aTrack.GetMomentumDirection());
  G4ThreeVector orth = momentum.orthogonal();
  G4ThreeVector newDir = momentum.rotate(orth, theta_charge);
  newDir.rotate(aTrack.GetMomentumDirection(), G4UniformRand()*2*pi);
  
  // FIXME:  Need to generate actual phonon!

  aParticleChange.ProposeMomentumDirection(newDir);
  aParticleChange.ProposeEnergy(T-qEnergy);
  aParticleChange.ProposeNonIonizingEnergyDeposit(qEnergy);
  ResetNumberOfInteractionLengthLeft();    
  
  return &aParticleChange;
}

G4bool G4CMPhLukeScattering::IsApplicable(const G4ParticleDefinition& aPD) {
  return(&aPD==G4CMPDriftHole::Definition() );
}

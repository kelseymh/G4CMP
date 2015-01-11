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
/// \file library/src/G4CMPhLukeScattering.cc
/// \brief Implementation of the G4CMPhLukeScattering class
//
// $Id$
//
// 20140325  Move time-step calculation to G4CMPProcessUtils
// 20140331  Add required process subtype code
// 20140509  Add run-time envvar to bias phonons
// 20141231  Rename "minimum step" function to ComputeMinTimeStep
// 20150111  Move envvar to G4CMPConfigManager

#include "G4CMPhLukeScattering.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftHole.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4PhononPolarization.hh"
#include "G4PhysicalConstants.hh"
#include "G4RandomDirection.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4VParticleChange.hh"
#include "Randomize.hh"
#include <iostream>
#include <fstream>


// Constructor and destructor

G4CMPhLukeScattering::G4CMPhLukeScattering(G4VProcess* stepper)
  : G4CMPVDriftProcess("hLukeScattering", fLukeScattering),
    stepLimiter(stepper) {
#ifdef G4CMP_DEBUG
  output.open("hLukePhononEnergies");
#endif
}

G4CMPhLukeScattering::~G4CMPhLukeScattering() {;}


// Physics

G4double 
G4CMPhLukeScattering::GetMeanFreePath(const G4Track& aTrack, G4double,
				      G4ForceCondition* condition) {
  *condition = Forced;		// In order to recompute MFP after TimeStepper
  
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
  G4StepPoint* postStepPoint = aStep.GetPostStepPoint();
  
#ifdef G4CMP_DEBUG
  G4cout << GetProcessName() << "::PostStepDoIt: Step limited by process "
	 << postStepPoint->GetProcessDefinedStep()->GetProcessName()
	 << G4endl;
#endif
  
  // Do nothing other than re-calculate mfp when step limit reached or leaving
  // volume
  if (postStepPoint->GetProcessDefinedStep()==stepLimiter
      || postStepPoint->GetStepStatus()==fGeomBoundary) {
    return &aParticleChange;
  }
  
  G4double velocity = postStepPoint->GetVelocity();
  G4double kmag = velocity*mc_h / hbar_Planck;

  G4double theta_phonon = MakePhononTheta(kmag, ksound_h);
  G4double phi_phonon = G4UniformRand()*twopi;
  G4double q = 2*(kmag*cos(theta_phonon)-ksound_h);

  G4double theta_charge = MakeRecoilTheta(kmag, ksound_h, theta_phonon);
  G4double phi_charge = fmod(pi+phi_phonon, twopi);

  // Rotate hole direction vector to recoil
  G4ThreeVector mdir = postStepPoint->GetMomentumDirection();
  G4ThreeVector newDir = mdir;
  newDir.rotate(mdir.orthogonal(), theta_charge);
  newDir.rotate(mdir, phi_charge);

  // Convert phonon pseudovector to real space
  G4ThreeVector qvec = q*mdir;
  qvec.rotate(mdir.orthogonal(), theta_phonon);
  qvec.rotate(mdir, phi_phonon);
  RotateToGlobalDirection(qvec);

  G4double Ephonon = MakePhononEnergy(kmag, ksound_h, theta_phonon);
#ifdef G4CMP_DEBUG
  output << Ephonon/eV << G4endl;
#endif

  // Create real phonon to be propagated, with random polarization
  static const G4double genLuke = G4CMPConfigManager::GetLukePhonons();
  if (genLuke > 0. && G4UniformRand() < genLuke) {
    G4Track* phonon = CreatePhonon(G4PhononPolarization::UNKNOWN,qvec,Ephonon);
    aParticleChange.SetNumberOfSecondaries(1);
    aParticleChange.AddSecondary(phonon);
  }

  G4double Etrack = postStepPoint->GetKineticEnergy();

  aParticleChange.ProposeMomentumDirection(newDir);
  aParticleChange.ProposeEnergy(Etrack-Ephonon);
  aParticleChange.ProposeNonIonizingEnergyDeposit(Ephonon);

  ResetNumberOfInteractionLengthLeft();
  return &aParticleChange;
}

G4bool G4CMPhLukeScattering::IsApplicable(const G4ParticleDefinition& aPD) {
  return (&aPD==G4CMPDriftHole::Definition());
}

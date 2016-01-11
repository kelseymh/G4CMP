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
/// \file library/src/G4CMPVLukeScattering.cc
/// \brief Implementation of the G4CMPVLukeScattering class
//
// $Id$
//
// 20150111  New base class for both electron and hole Luke processes
// 20150122  Use verboseLevel instead of compiler flag for debugging

#include "G4CMPLukeScattering.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4CMPTrackInformation.hh"
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

G4CMPLukeScattering::G4CMPLukeScattering(G4VProcess* stepper)
  : G4CMPVDriftProcess("G4CMPLukeScattering", fLukeScattering),
    stepLimiter(stepper) {
#ifdef G4CMP_DEBUG
  output.open("LukePhononEnergies");
#endif
}

G4CMPLukeScattering::~G4CMPLukeScattering() {
#ifdef G4CMP_DEBUG
  output.close();
#endif
}


// Physics

G4double 
G4CMPLukeScattering::GetMeanFreePath(const G4Track& aTrack, G4double,
                                     G4ForceCondition* condition) {
  *condition = Forced;		// In order to recompute MFP after TimeStepper

  G4CMPTrackInformation* trackInfo = static_cast<G4CMPTrackInformation*>(
    aTrack.GetAuxiliaryTrackInformation(fPhysicsModelID));

  G4double velocity = GetVelocity(aTrack);
  G4double kmag = GetLocalWaveVector(aTrack).mag();
  G4double l0 = trackInfo->GetScatterLength();
  G4double kSound = CalculateKSound(trackInfo);

  if (verboseLevel > 1) {
    G4cout << "LukeScattering v = " << velocity/m*s << " kmag = " << kmag*m
	   << G4endl;
  }

  if (kmag <= kSound) return DBL_MAX;
  
  // Time step corresponding to Mach number (avg. time between radiations)
  G4double dtau = ChargeCarrierTimeStep(kmag/kSound, l0);
  G4double mfp = dtau * velocity;

  if (verboseLevel > 1) {
    G4cout << "LukeScattering MFP = " << mfp << G4endl;
  }

  return mfp;
}


G4VParticleChange* G4CMPLukeScattering::PostStepDoIt(const G4Track& aTrack,
                                                     const G4Step& aStep) {
  aParticleChange.Initialize(aTrack); 
  G4StepPoint* postStepPoint = aStep.GetPostStepPoint();
  
  // Do nothing other than re-calculate mfp when step limit reached or leaving
  // volume
  if (verboseLevel > 1) {
    G4cout << GetProcessName() << "::PostStepDoIt: Step limited by process "
	   << postStepPoint->GetProcessDefinedStep()->GetProcessName()
	   << G4endl;
  }

  // Do nothing other than re-calculate mfp when step limit reached or leaving volume
  if (postStepPoint->GetStepStatus()==fGeomBoundary ||
      postStepPoint->GetProcessDefinedStep()==stepLimiter) {
    return &aParticleChange;
  }

  G4CMPTrackInformation* trackInfo = static_cast<G4CMPTrackInformation*>(
    aTrack.GetAuxiliaryTrackInformation(fPhysicsModelID));

  G4ThreeVector ktrk = GetLocalWaveVector(aTrack);
  G4double kmag = ktrk.mag();
  G4double kSound = CalculateKSound(trackInfo);

  // Sanity check: this should have been done in MFP already
  if (kmag <= kSound) return &aParticleChange;

  if (verboseLevel > 1) {
    G4cout << "p (post-step) = " << postStepPoint->GetMomentum()
	   << "\np_mag = " << postStepPoint->GetMomentum().mag()
	   << "\nktrk = " << ktrk
     << "\nkmag = " << kmag << " k/ks = " << kmag/kSound
     << "\nacos(ks/k) = " << acos(kSound/kmag) << G4endl;
  }

  G4double theta_phonon = MakePhononTheta(kmag, kSound);
  G4double phi_phonon   = G4UniformRand()*twopi;
  G4double q = 2*(kmag*cos(theta_phonon)-kSound);

  // Sanity check for phonon production: should be forward, like Cherenkov
  if (theta_phonon>acos(kSound/kmag) || theta_phonon>halfpi) {
    G4cerr << GetProcessName() << " ERROR: Phonon production theta_phonon "
     << theta_phonon << " exceeds cone angle " << acos(kSound/kmag)
	   << G4endl;
    return &aParticleChange;
  }
  
  // Generate phonon momentum vector
  G4ThreeVector kdir = ktrk.unit();
  G4ThreeVector qvec = q*kdir;
  qvec.rotate(kdir.orthogonal(), theta_phonon);
  qvec.rotate(kdir, phi_phonon);

  G4double Ephonon = MakePhononEnergy(kmag, kSound, theta_phonon);
#ifdef G4CMP_DEBUG
  output << Ephonon/eV << G4endl;
#endif

  // Get recoil wavevector, convert to new momentum
  G4ThreeVector k_recoil = ktrk - qvec;

  if (verboseLevel > 1) {
    G4cout << "theta_phonon = " << theta_phonon
	   << " phi_phonon = " << phi_phonon
	   << "\nq = " << q << "\nqvec = " << qvec << "\nEphonon = " << Ephonon
	   << "\nk_recoil = " << k_recoil
	   << "\nk_recoil-mag = " << k_recoil.mag()
	   << G4endl;
  }

  // Create real phonon to be propagated, with random polarization
  static const G4double genLuke = G4CMPConfigManager::GetGenPhonons();
  //if (genLuke > 0. && G4UniformRand() < genLuke) {
  if (false) {
    MakeGlobalPhononK(qvec);  		// Convert phonon vector to real space

    G4Track* phonon = CreatePhonon(G4PhononPolarization::UNKNOWN,qvec,Ephonon);
    aParticleChange.SetNumberOfSecondaries(1);
    aParticleChange.AddSecondary(phonon);
  }

  MakeGlobalRecoil(k_recoil);		// Converts wavevector to momentum
  FillParticleChange(GetValleyIndex(aTrack), k_recoil);

  aParticleChange.ProposeNonIonizingEnergyDeposit(Ephonon);
  ResetNumberOfInteractionLengthLeft();
  return &aParticleChange;
}

G4double
G4CMPLukeScattering::CalculateKSound(const G4CMPTrackInformation* trackInfo) {
  return theLattice->GetSoundSpeed()*trackInfo->GetEffectiveMass()/hbar_Planck;
}

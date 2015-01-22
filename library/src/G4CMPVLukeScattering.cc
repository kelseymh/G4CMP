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

#include "G4CMPVLukeScattering.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftElectron.hh"
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

G4CMPVLukeScattering::G4CMPVLukeScattering(const G4String& name,
					   const G4ParticleDefinition* carrier,
					   G4VProcess* stepper)
  : G4CMPVDriftProcess(name+"Scattering", fLukeScattering),
    shortName(name), stepLimiter(stepper), theCarrier(carrier),
    theKsound(ksound_e), theL0(l0_e) {
#ifdef G4CMP_DEBUG
  output.open((shortName+"PhononEnergies").c_str());
#endif
}

G4CMPVLukeScattering::~G4CMPVLukeScattering() {
#ifdef G4CMP_DEBUG
  output.close();
#endif
}


// Configuration

void G4CMPVLukeScattering::LoadDataForTrack(const G4Track* track) {
  if (track->GetDefinition() != theCarrier) {
    G4cerr << GetProcessName() << " ERROR:  Track type "
	   << track->GetDefinition()->GetParticleName() << " not valid" << G4endl;
    return;
  }

  G4CMPVDriftProcess::LoadDataForTrack(track);

  if (theCarrier == G4CMPDriftHole::Definition()) {
    theKsound = ksound_h;
    theL0 = l0_h;
  } else if (theCarrier == G4CMPDriftElectron::Definition()) {
    theKsound = ksound_e;
    theL0 = l0_e;
  }
}

G4bool G4CMPVLukeScattering::IsApplicable(const G4ParticleDefinition& aPD) {
  return (&aPD==theCarrier);
}

G4double G4CMPVLukeScattering::GetWaveNumber(const G4Track& aTrack) const {
  return GetWaveVector(aTrack).mag();
}


// Physics

G4double 
G4CMPVLukeScattering::GetMeanFreePath(const G4Track& aTrack, G4double,
				      G4ForceCondition* condition) {
  *condition = Forced;		// In order to recompute MFP after TimeStepper

  G4double velocity = GetVelocity(aTrack);
  G4double kmag = GetWaveNumber(aTrack);
  
  if (verboseLevel > 1) {
    G4cout << shortName << " v = " << velocity/m*s << " kmag = " << kmag*m
	   << G4endl;
  }

  if (kmag<=theKsound) return DBL_MAX;
  
  // Time step corresponding to Mach number (avg. time between radiations)
  G4double dtau = ChargeCarrierTimeStep(kmag/theKsound, theL0);
  G4double mfp = dtau * velocity;

  if (verboseLevel > 1) {
    G4cout << shortName << " MFP = " <<  mfp <<  G4endl;
  }

  return mfp;
}


G4VParticleChange* G4CMPVLukeScattering::PostStepDoIt(const G4Track& aTrack,
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

  G4double kmag = GetWaveNumber(aTrack);
  G4ThreeVector ktrk = GetWaveVector(aTrack);

  // Sanity check: this should have been done in MFP already
  if (kmag <= theKsound) return &aParticleChange;

  if (verboseLevel > 1) {
    G4cout << "p (post-step) = " << postStepPoint->GetMomentum()
	   << "\np_mag = " << postStepPoint->GetMomentum().mag()
	   << "\nktrk = " << ktrk
	   << "\nkmag = " << kmag << " k/ks = " << kmag/theKsound
	   << "\nacos(ks/k) = " << acos(theKsound/kmag) << G4endl;
  }

  G4double theta_phonon = MakePhononTheta(kmag, theKsound);
  G4double phi_phonon   = G4UniformRand()*twopi;
  G4double q = 2*(kmag*cos(theta_phonon)-theKsound);

  // Sanity check for phonon production: should be forward, like Cherenkov
  if (theta_phonon>acos(theKsound/kmag) || theta_phonon>halfpi) {
    G4cerr << GetProcessName() << " ERROR: Phonon production theta_phonon "
	   << theta_phonon << " exceeds cone angle " << acos(theKsound/kmag)
	   << G4endl;
    return &aParticleChange;
  }
  
  // Generate phonon momentum vector
  G4ThreeVector kdir = ktrk.unit();
  G4ThreeVector qvec = q*kdir;
  qvec.rotate(kdir.orthogonal(), theta_phonon);
  qvec.rotate(kdir, phi_phonon);

  G4double Ephonon = MakePhononEnergy(kmag, theKsound, theta_phonon);
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
  static const G4double genLuke = G4CMPConfigManager::GetLukePhonons();
  if (genLuke > 0. && G4UniformRand() < genLuke) {
    MakeGlobalPhonon(qvec);  		// Convert phonon vector to real space

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

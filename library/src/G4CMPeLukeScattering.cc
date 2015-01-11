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
/// \file library/src/G4CMPeLukeScattering.cc
/// \brief Implementation of the G4CMPeLukeScattering class
//
// $Id$
//
// 20140325  Move time-step calculation to G4CMPProcessUtils
// 20140331  Add required process subtype code
// 20140415  Add run-time flag to select valley vs. H-V kinematics
// 20140430  Compute kinematics using mass tensor; prepare to create phonons
// 20140509  Remove valley vs. H-V flag; add run-time envvar to bias phonons
// 20140521  Remove momentum-check loop; energy conservation is enforced
// 20140903  Get Etrack using valley kinematics, _not_ track or stepPoint
// 201411??  R.Agnese -- Merge functionality from TimeStepper here
// 20141231  Rename "minimum step" function to ComputeMinTimeStep
// 20150106  Move envvar to G4CMPConfigManager

#include "G4CMPeLukeScattering.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPValleyTrackMap.hh"
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

G4CMPeLukeScattering::G4CMPeLukeScattering(G4VProcess* sLim)
  : G4CMPVDriftProcess("eLukeScattering", fLukeScattering),
    stepLimiter(sLim) {
#ifdef G4CMP_DEBUG
  output.open("eLukePhononEnergies");
#endif
}

G4CMPeLukeScattering::~G4CMPeLukeScattering() {;}


// Physics

G4double G4CMPeLukeScattering::GetMeanFreePath(const G4Track& aTrack,
					       G4double,
					       G4ForceCondition* condition) {
  *condition = Forced;

  G4int iv = GetValleyIndex(aTrack);
  G4ThreeVector p_local = GetLocalMomentum(aTrack);
  G4ThreeVector v = theLattice->MapPtoV_el(iv, p_local);
  G4ThreeVector k_HV = theLattice->MapPtoK_HV(iv, p_local);
  G4ThreeVector k_valley = theLattice->MapPtoK_valley(iv, p_local);

  G4double kmag = k_HV.mag();

#ifdef G4CMP_DEBUG
  G4cout << "eLuke v = " << v.mag()/m*s << " kmag = " << kmag*m
	 << "\v = " << v << "\nk_HV = " << k_HV
	 << "\nk_valley = " << k_valley << G4endl;
#endif

  if (kmag<=ksound_e) return DBL_MAX;
 
  // Time step corresponding to Mach number (avg. time between radiations)
  G4double dtau = ChargeCarrierTimeStep(kmag/ksound_e, l0_e);
  
  G4double mfp = dtau * v.mag();
#ifdef G4CMP_DEBUG
  G4cout << "eLuke MFP = " << mfp/m << G4endl;
#endif

  return mfp;
}

G4VParticleChange* G4CMPeLukeScattering::PostStepDoIt(const G4Track& aTrack,
						      const G4Step& aStep) {
  aParticleChange.Initialize(aTrack); 
  G4StepPoint* postStepPoint = aStep.GetPostStepPoint();
  
  // Do nothing other than re-calculate mfp when step limit reached or leaving
  // volume
#ifdef G4CMP_DEBUG
  G4cout << GetProcessName() << "::PostStepDoIt: Step limited by process "
	 << postStepPoint->GetProcessDefinedStep()->GetProcessName()
	 << G4endl;
#endif

  // Do nothing other than re-calculate mfp when step limit reached or leaving volume
  if (postStepPoint->GetStepStatus()==fGeomBoundary ||
      postStepPoint->GetProcessDefinedStep()==stepLimiter) {
    return &aParticleChange;
  }

  G4int iv = GetValleyIndex(aTrack);
  G4ThreeVector p = GetLocalDirection(postStepPoint->GetMomentum());
  G4ThreeVector k_valley = theLattice->MapPtoK_valley(iv, p);
  G4ThreeVector k_HV = theLattice->MapPtoK_HV(iv, p);

  G4double kmag = k_HV.mag();

  if (kmag <= ksound_e) return &aParticleChange;

#ifdef G4CMP_DEBUG
  G4cout << "p (post-step) = " << p << "\np_mag = " << p.mag()
	 << "\nk_HV = " << k_HV << " k_valley = " << k_valley
	 << "\nkmag = " << kmag << " k/ks = " << kmag/ksound_e
	 << "\nacos(ks/k)   = " << acos(ksound_e/kmag) << G4endl;
#endif

  // Compute kinematics in Herring-Vogt space, where electron mass is scalar
  G4double theta_phonon=0., phi_phonon=0., q=0., Enew=0.;
  G4ThreeVector kdir, qvec, k_recoil, p_new;

  G4double Etrack = theLattice->MapPtoEkin(iv, p);

  // Iterate phonon generation to avoid energy violations
  const G4int maxTries = 1000;
  G4int nTries = 0;
  do {
    nTries++;

    // Construct phonon pseudovector in Herring-Vogt space
    theta_phonon = MakePhononTheta(kmag, ksound_e);
    phi_phonon   = G4UniformRand()*twopi;
    q = 2*(kmag*cos(theta_phonon)-ksound_e);
    
    // Sanity check for phonon production: should be forward, like Cherenkov
    if (theta_phonon>acos(ksound_e/kmag) || theta_phonon>halfpi) {
      G4cerr << "ERROR: Phonon production theta_phonon " << theta_phonon
	     << " exceeds cone angle " << acos(ksound_e/kmag) << G4endl;
      continue;
    }
    
    kdir = k_HV.unit();			// Get phonon and recoil vectors
    qvec = q*kdir;
    qvec.rotate(kdir.orthogonal(), theta_phonon);
    qvec.rotate(kdir, phi_phonon);

    // Get recoil wavevector in HV space, convert to new momentum
    k_recoil = k_HV - qvec;
    p_new = theLattice->MapK_HVtoP(iv, k_recoil);

    Enew = theLattice->MapPtoEkin(iv, p_new);

    // Sanity check: electron should have lost energy in recoil
  } while (Enew >= Etrack && nTries < maxTries);

#ifdef G4CMP_DEBUG
  G4cout << GetProcessName() << " Phonon generation " << nTries << " attempts"
	 << G4endl;
#endif

  // Sanity check: electron should have lost energy in recoil
  if (Enew > Etrack) {
    G4cerr << GetProcessName() << " ERROR: Recoil energy exceeds input after "
	     << nTries << " attempts." << G4endl;
    return &aParticleChange;
  }

  // Convert phonon pseudovector to real space
  G4double Ephonon = MakePhononEnergy(kmag, ksound_e, theta_phonon);
#ifdef G4CMP_DEBUG
  output << Ephonon/eV << G4endl;
#endif

  qvec = theLattice->MapK_HVtoK_valley(iv, qvec);
  RotateToGlobalDirection(qvec);

  // Create real phonon to be propagated, with random polarization
  static const G4double genLuke = G4CMPConfigManager::GetLukePhonons();
  if (genLuke > 0. && G4UniformRand() < genLuke) {
    G4Track* phonon = CreatePhonon(G4PhononPolarization::UNKNOWN,qvec,Ephonon);
    aParticleChange.SetNumberOfSecondaries(1);
    aParticleChange.AddSecondary(phonon);
  }

#ifdef G4CMP_DEBUG
  G4cout << "\ntheta_phonon = " << theta_phonon
	 << " phi_phonon = " << phi_phonon
	 << "\nq = " << q << "\nqvec = " << qvec << "\nEphonon = " << Ephonon
	 << "\nk_recoil = " << k_recoil << "\nk_recoil-mag = " << k_recoil.mag()
	 << "\np_new    = " << p_new << "\np_new_mag = " << p_new.mag()
	 << "\nEtrack = " << Etrack << " Enew = " << Enew
	 << G4endl;
#endif

  // Compute energy and "effective mass" of recoiling electron
  RotateToGlobalDirection(p_new);	// Put into global frame for tracks
  SetNewKinematics(iv, p_new);

  aParticleChange.ProposeNonIonizingEnergyDeposit(Ephonon);

  ResetNumberOfInteractionLengthLeft();
  return &aParticleChange;
}

G4bool G4CMPeLukeScattering::IsApplicable(const G4ParticleDefinition& aPD) {
  return (&aPD==G4CMPDriftElectron::Definition());
}

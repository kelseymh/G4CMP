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
// 20140325  Move time-step calculation to G4CMPProcessUtils
// 20140331  Add required process subtype code
// 20140415  Add run-time flag to select valley vs. H-V kinematics
// 20140430  Compute kinematics using mass tensor; prepare to create phonons

#include "G4CMPeLukeScattering.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPValleyTrackMap.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhononPolarization.hh"
#include "G4RandomDirection.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4VParticleChange.hh"
#include "Randomize.hh"
#include <fstream>
#include <iostream>


G4CMPeLukeScattering::G4CMPeLukeScattering(G4VProcess* sLim)
  : G4CMPVDriftProcess("eLukeScattering", fLukeScattering),
    stepLimiter(sLim) {;}

G4CMPeLukeScattering::~G4CMPeLukeScattering() {;}


G4double G4CMPeLukeScattering::GetMeanFreePath(const G4Track& aTrack,
					       G4double,
					       G4ForceCondition* condition) {
  *condition = Forced;

  G4int iv = GetValleyIndex(aTrack);
  G4ThreeVector v = theLattice->MapPtoV_el(iv,GetLocalMomentum(aTrack));
  G4ThreeVector k_HV = theLattice->MapPtoK_HV(iv,GetLocalMomentum(aTrack));
  G4double kmag = k_HV.mag();

#ifdef G4CMP_DEBUG
  G4cout << "eLuke v = " << v.mag()/m*s << " kmag = " << kmag*m
	 << "\nv = " << v << "\nk_HV = " << k_HV
	 << "\nk_valley = "
	 << theLattice->MapPtoK_valley(iv,GetLocalMomentum(aTrack))
	 << G4endl;
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

  if ( (postStepPoint->GetProcessDefinedStep()==stepLimiter)
       || (postStepPoint->GetStepStatus()==fGeomBoundary) ) {
    return &aParticleChange;
  }

  G4int iv = GetValleyIndex(aTrack);
  G4ThreeVector p = GetLocalDirection(postStepPoint->GetMomentum());
  G4ThreeVector k_valley = theLattice->MapPtoK_valley(iv, p);
  G4ThreeVector k_HV = theLattice->MapPtoK_HV(iv, p);
  G4double kmag = k_HV.mag();

#ifdef G4CMP_DEBUG
  G4cout << "p (post-step) = " << p << "\np_mag = " << p.mag()
	 << "\nk_HV = " << k_HV
	 << "\nkmag = " << kmag << " k/ks = " << kmag/ksound_e
	 << "\nacos(ks/k)   = " << acos(ksound_e/kmag) << G4endl;
#endif

  // Compute kinematics in Herring-Vogt space, where electron mass is scalar
  G4double theta_phonon=0., phi_phonon=0., q=0.;
  G4ThreeVector kdir, qvec, k_recoil, p_new;

  // Iterate until final momentum (magnitude) less than track
  G4int ntries = 0;
  do {
    ntries++;
#ifdef G4CMP_DEBUG
    if (ntries%10000 == 0) G4cout << "... generator trial " << ntries << G4endl;
#endif
    
    // Construct phonon pseudovector in Herring-Vogt space
    theta_phonon = MakePhononTheta(kmag, ksound_e);
    phi_phonon   = G4UniformRand()*twopi;
    q = 2*(kmag*cos(theta_phonon)-ksound_e);

    // Sanity check for phonon production: should be forward, like Cherenkov
    if (theta_phonon>acos(ksound_e/kmag) || theta_phonon>halfpi) {
      G4cerr << "ERROR: Phonon production theta_phonon " << theta_phonon
	     << " exceeds cone angle " << acos(ksound_e/kmag) << G4endl;
    }
    
    kdir = k_HV.unit();			// Get phonon and recoil vectors
    qvec = q*kdir;
    qvec.rotate(kdir.orthogonal(), theta_phonon);
    qvec.rotate(kdir, phi_phonon);
    k_recoil = k_HV - qvec;
    p_new = theLattice->MapK_HVtoP(iv, k_recoil);
  } while ((p_new.mag()-p.mag())/p.mag() > -1e-7);	// Check conservation

#ifdef G4CMP_DEBUG
  G4cout << "Phonon generation required " << ntries << " trials" << G4endl;
#endif

  // Convert phonon pseudovector to real space
  G4double Ephonon = MakePhononEnergy(kmag, ksound_e, theta_phonon);
  qvec = theLattice->MapK_HVtoK_valley(iv, qvec);
  RotateToGlobalDirection(qvec);

  /*****
  // Create real phonon to be propagated
  G4Track* phonon = CreatePhonon(G4PhononPolarization::Long, qvec, Ephonon);
  aParticleChange.SetNumberOfSecondaries(1);
  aParticleChange.AddSecondary(phonon);
  *****/

#ifdef G4CMP_DEBUG
  G4double Etrack = theLattice->MapPtoEkin(iv, p);
  G4double Enew = theLattice->MapPtoEkin(iv, p_new);

  G4cout << "\ntheta_phonon = " << theta_phonon
	 << " phi_phonon = " << phi_phonon
	 << "\nq = " << q << "\nqvec = " << qvec << "\nEphonon = " << Ephonon
	 << "\nk_recoil = " << k_recoil << "\nk_recoil-mag = " << k_recoil.mag()
	 << "\np_new     = " << p_new << "\np_new_mag = " << p_new.mag()
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

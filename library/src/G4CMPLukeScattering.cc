/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPVLukeScattering.cc
/// \brief Implementation of the G4CMPVLukeScattering class
//
// $Id$
//
// 20150111  New base class for both electron and hole Luke processes
// 20150122  Use verboseLevel instead of compiler flag for debugging
// 20160624  Use GetTrackInfo() accessor
// 20160830  Replace direct use of G4CMP_MAKE_PHONONS with ChooseWeight
// 20161114  Use new DriftTrackInfo
// 20170602  Use G4CMPUtils for track identity functions
// 20170802  Use G4CMP_LUKE_SAMPLE biasing with ChooseWeight()

#include "G4CMPLukeScattering.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4CMPDriftTrackInfo.hh"
#include "G4CMPSecondaryUtils.hh"
#include "G4CMPTrackUtils.hh"
#include "G4CMPUtils.hh"
#include "G4ExceptionSeverity.hh"
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
  if (!output.good()) {
    G4Exception("G4LatticeReader::MakeLattice", "Lattice001",
		FatalException, "Unable to open LukePhononEnergies");
  }
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

  auto trackInfo = G4CMP::GetTrackInfo<G4CMPDriftTrackInfo>(aTrack);
  const G4LatticePhysical* lat = trackInfo->Lattice();
  G4double kmag = 0.; G4double l0 = 0.; G4double mass = 0.;
  if (G4CMP::IsElectron(aTrack)) {
    kmag = lat->MapV_elToK_HV(GetValleyIndex(aTrack),
                              GetLocalVelocityVector(aTrack)).mag();
    l0 = lat->GetElectronScatter();
    mass = lat->GetElectronMass();
  } else if (G4CMP::IsHole(aTrack)) {
    kmag = GetLocalWaveVector(aTrack).mag();
    l0 = lat->GetHoleScatter();
    mass = lat->GetHoleMass();
  } else {
    G4Exception("G4CMPLukeScattering::GetMeanFreePath", "Luke001",
                EventMustBeAborted, "Unknown charge carrier");
  }

  G4double kSound = lat->GetSoundSpeed() * mass / hbar_Planck;
  G4double velocity = GetVelocity(aTrack);

  if (verboseLevel > 1) {
    G4cout << "LukeScattering v = " << velocity/m*s << " kmag = " << kmag*m
	   << G4endl;
  }

  if (kmag <= kSound) return DBL_MAX;
  
  // Time step corresponding to Mach number (avg. time between radiations)
  G4double dtau = ChargeCarrierTimeStep(kmag/kSound, l0);
  G4double mfp = dtau * velocity;

  if (verboseLevel > 1) {
    G4cout << "LukeScattering Time Step = " << dtau << "\n"
           << "LukeScattering MFP = " << mfp << G4endl;
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

  auto trackInfo = G4CMP::GetTrackInfo<G4CMPDriftTrackInfo>(aTrack);
  const G4LatticePhysical* lat = trackInfo->Lattice();

  G4ThreeVector ktrk(0.);
  G4double mass = 0.;
  if (IsElectron()) {
    ktrk = lat->MapV_elToK_HV(GetValleyIndex(aTrack),
                              GetLocalVelocityVector(aTrack));
    mass = lat->GetElectronMass();
  } else if (IsHole()) {
    ktrk = GetLocalWaveVector(aTrack);
    mass = lat->GetHoleMass();
  } else {
    G4Exception("G4CMPLukeScattering::PostStepDoIt", "Luke002",
                EventMustBeAborted, "Unknown charge carrier");
  }

  G4double kmag = ktrk.mag();
  G4double kSound = lat->GetSoundSpeed() * mass / hbar_Planck;

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
  // If phonon is not created, register the energy as deposited
  G4double weight =
    G4CMP::ChoosePhononWeight(G4CMPConfigManager::GetLukeSampling());

  if (weight > 0.) {
    MakeGlobalPhononK(qvec);  		// Convert phonon vector to real space

    G4Track* phonon = G4CMP::CreatePhonon(aTrack.GetTouchable(),
                                          G4PhononPolarization::UNKNOWN,
                                          qvec,Ephonon,
                                          aTrack.GetGlobalTime(),
                                          aTrack.GetPosition());
    // Secondary's weight has to be multiplicative with its parent's
    phonon->SetWeight(aTrack.GetWeight() * weight);

    aParticleChange.SetSecondaryWeightByProcess(true);
    aParticleChange.SetNumberOfSecondaries(1);
    aParticleChange.AddSecondary(phonon);
  } else {
    aParticleChange.ProposeNonIonizingEnergyDeposit(Ephonon);
  }

  MakeGlobalRecoil(k_recoil);		// Converts wavevector to momentum
  FillParticleChange(GetValleyIndex(aTrack), k_recoil);
  ResetNumberOfInteractionLengthLeft();
  return &aParticleChange;
}

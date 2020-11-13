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
// 20170805  Use scattering-rate model
// 20170907  Make process non-forced; check only for boundary crossing
// 20170928  Hide "output" usage behind verbosity check, as well as G4CMP_DEBUG
// 20180827  Add debugging output with weight calculation.
// 20190816  Add flag to track secondary phonons immediately (c.f. G4Cerenkov)
// 20201109  Modify debugging output file with additional information, move
//		debugging output creation to PostStepDoIt to allows settting
//		process verbosity via macro commands.

#include "G4CMPLukeScattering.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4CMPDriftTrackInfo.hh"
#include "G4CMPLukeEmissionRate.hh"
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
    stepLimiter(stepper), secondariesFirst(true) {
  UseRateModel(new G4CMPLukeEmissionRate);
}

G4CMPLukeScattering::~G4CMPLukeScattering() {
#ifdef G4CMP_DEBUG
  if (output.is_open()) output.close();
#endif
}


// Physics

G4VParticleChange* G4CMPLukeScattering::PostStepDoIt(const G4Track& aTrack,
                                                     const G4Step& aStep) {
  aParticleChange.Initialize(aTrack); 
  G4StepPoint* postStepPoint = aStep.GetPostStepPoint();
  
  if (verboseLevel > 1) {
    G4cout << GetProcessName() << "::PostStepDoIt: Step limited by process "
           << postStepPoint->GetProcessDefinedStep()->GetProcessName()
           << G4endl;
  }

  // Don't do anything at a volume boundary
  if (postStepPoint->GetStepStatus()==fGeomBoundary) {
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }

#ifdef G4CMP_DEBUG
  if (verboseLevel && !output.is_open()) {
    output.open("LukePhononEnergies");
    if (!output.good()) {
      G4Exception("G4LatticeReader::MakeLattice", "Lattice001",
		  FatalException, "Unable to open LukePhononEnergies");
    }

    output << "Track Type,Track Energy [eV],WaveVector,Phonon Theta,"
	   << "Phonon Energy [eV],Recoil WaveVector,Final Energy [eV]"
	   << std::endl;
  }
#endif

  // Collect ancillary information needed for kinematics
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
    return &aParticleChange;
  }

  G4double kmag = ktrk.mag();
  G4double kSound = lat->GetSoundSpeed() * mass / hbar_Planck;

  // Sanity check: this should have been done in MFP already
  if (kmag <= kSound) return &aParticleChange;

  if (verboseLevel > 1) {
    G4cout << "p (post-step) = " << postStepPoint->GetMomentum()
	   << "\np_mag = " << postStepPoint->GetMomentum().mag()
	   << "\nktrk = " << ktrk << " kmag = " << kmag
	   << "\nk/ks = " << kmag/kSound
	   << " acos(ks/k) = " << acos(kSound/kmag) << G4endl;
  }

  G4double theta_phonon = MakePhononTheta(kmag, kSound);
  G4double phi_phonon   = G4UniformRand()*twopi;
  G4double q = 2*(kmag*cos(theta_phonon)-kSound);

  if (verboseLevel > 1) {
    G4cout << "theta_phonon = " << theta_phonon
           << " phi_phonon = " << phi_phonon << " q = " << q << G4endl;
  }

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

  if (verboseLevel > 1) {
    G4cout << "qvec = " << qvec
	   << "\nktrk.qvec = " << ktrk.dot(qvec)/(kmag*q)
	   << " ktr-qvec angle " << acos(ktrk.dot(qvec)/(kmag*q))
           << G4endl;
  }

  // Get recoil wavevector, convert to new momentum
  G4ThreeVector k_recoil = ktrk - qvec;

  // Phonon vector and energy should be in local frame, not H-V
  if (IsElectron()) 
    qvec = lat->MapK_HVtoK(GetValleyIndex(aTrack), qvec);

  G4double Ephonon = MakePhononEnergy(qvec.mag());

  // Sanity check for phonon production: can't exceed charge's energy
  if (Ephonon >= GetKineticEnergy(aTrack)) {
    G4cerr << GetProcessName() << "ERROR: Phonon production Ephonon "
	   << Ephonon/eV << " eV exceeds charge carrier energy "
	   << GetKineticEnergy(aTrack)/eV << " eV" << G4endl;
    return &aParticleChange;
  }

  if (verboseLevel > 1) {
    G4cout << "q(HV) = " << q << " q(local) = " << qvec.mag()
	   << "\nEphonon = " << Ephonon
           << "\nk_recoil = " << k_recoil
           << " k_recoil-mag = " << k_recoil.mag()
           << G4endl;
  }

#ifdef G4CMP_DEBUG
  if (output.good()) {
    output << aTrack.GetDefinition()->GetParticleName() << ","
	   << GetKineticEnergy(aTrack)/eV << "," << kmag << ","
	   << theta_phonon << "," << Ephonon/eV << "," << k_recoil.mag()
	   << ",";
  }
#endif

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
    if (verboseLevel>1) {
      G4cout << "phonon wt " << phonon->GetWeight()
	     << " : track " << aTrack.GetTrackID()
	     << " wt " << aTrack.GetWeight()
	     << "  thrown wt " << weight << G4endl;
    }

    aParticleChange.SetSecondaryWeightByProcess(true);
    aParticleChange.SetNumberOfSecondaries(1);
    aParticleChange.AddSecondary(phonon);

    // If user wants to track phonons immediately, put track back on stack
    if (secondariesFirst && aTrack.GetTrackStatus() == fAlive)
      aParticleChange.ProposeTrackStatus(fSuspend);
  } else {
    aParticleChange.ProposeNonIonizingEnergyDeposit(Ephonon);
  }

  MakeGlobalRecoil(k_recoil);		// Converts wavevector to momentum
  FillParticleChange(GetValleyIndex(aTrack), k_recoil);

#ifdef G4CMP_DEBUG
  if (output.good()) output << aParticleChange.GetEnergy()/eV << std::endl;
#endif

  ClearNumberOfInteractionLengthLeft();
  return &aParticleChange;
}

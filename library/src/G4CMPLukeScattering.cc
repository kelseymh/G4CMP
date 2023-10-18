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
// 20201119  Put kinematics selection in accept/reject loop insted of quitting.
// 20201124  Include track momenta (before and after) in diagnostic output.
// 20201207  For electrons, add energy conservation check in accept/reject.
// 20210225  Improve kinematics: do not transform the phonon vector.  Leave
//		in place code for electrons to assign final-state to valley
//		closest to momentum direction.  Commented out now, as it leads
//		to non-physical reduction of total Luke emission.
// 20220907  G4CMP-316 -- Pass track into CreatePhonon instead of touchable.

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

    output << "Track ID, Track Type,Track Weight,Track Energy [eV],Track Momentum [eV],WaveVector,"
	   << "Phonon Theta,Phonon Energy [eV],Phonon Weight,Recoil WaveVector,"
	   << "Final Energy [eV],Final Momentum [eV]"
	   << std::endl;
  }
#endif

  // For convenience in diagnostic output below
  const G4String& trkName = aTrack.GetDefinition()->GetParticleName();

  // Collect ancillary information needed for kinematics
  auto trackInfo = G4CMP::GetTrackInfo<G4CMPDriftTrackInfo>(aTrack);
  const G4LatticePhysical* lat = trackInfo->Lattice();

  G4int iValley = GetValleyIndex(aTrack);	// Doesn't change valley

  // NOTE: Track kinematics include post-step acceleration from E-field
  G4ThreeVector ptrk = GetLocalMomentum(aTrack);
  G4ThreeVector ktrk(0.);
  G4double mass = 0.;
  if (IsElectron()) {
    ktrk = lat->MapV_elToK_HV(iValley, GetLocalVelocityVector(aTrack));
    mass = lat->GetElectronMass();
  } else if (IsHole()) {
    ktrk = GetLocalWaveVector(aTrack);
    mass = lat->GetHoleMass();
  } else {
    G4Exception("G4CMPLukeScattering::PostStepDoIt", "Luke002",
                EventMustBeAborted, "Unknown charge carrier");
    return &aParticleChange;
  }

  G4ThreeVector kdir = ktrk.unit();
  G4double kmag = ktrk.mag();
  G4double kSound = lat->GetSoundSpeed() * mass / hbar_Planck;

  // Sanity check: this should have been done in MFP already
  if (kmag <= kSound) return &aParticleChange;

  if (verboseLevel > 1) {
    G4cout << "p_global = " << GetGlobalMomentum(aTrack)
	   << " p_mag " << GetGlobalMomentum(aTrack).mag() << G4endl
	   << " p_local = " << ptrk << " p_mag " << ptrk.mag() << G4endl;

    if (IsElectron()) {
      G4ThreeVector kvalley = theLattice->MapPtoK_valley(iValley, ptrk);
      G4ThreeVector pvalley = kvalley * hbarc;
      G4cout << " valley " << iValley << " along "
	     << theLattice->GetValleyAxis(iValley) << G4endl
	     << " p_valley = " << pvalley << " p_mag " << pvalley.mag()
	     << G4endl
	     << " k_valley = " << kvalley << " k_mag " << kvalley.mag()
	     << G4endl;
    }

    G4cout << " ktrk = " << ktrk << " kmag " << kmag << G4endl
	   << " k/ks = " << kmag/kSound << " acos(ks/k) = " << acos(kSound/kmag)
	   << G4endl;
  }

  // Final state kinematics, generated in accept/reject loop below
  G4double theta_phonon=0, phi_phonon=0, q=0, Ephonon=0, Erecoil=0;
  G4ThreeVector qvec, k_recoil, precoil;	// Outgoing wave vectors
  G4int newValley = iValley;			// Test for valley change

  // Iterate to avoid non-physical phonon emission
  const G4int maxThrows = 1000;		// Avoids potential infinite loop
  G4bool goodThrow = false;
  G4int iThrow = 0;
  while (!goodThrow && iThrow++ < maxThrows) {
    theta_phonon = MakePhononTheta(kmag, kSound);
    phi_phonon   = G4UniformRand()*twopi;
    q = 2*(kmag*cos(theta_phonon)-kSound);

    if (verboseLevel > 1) {
      G4cout << " theta_phonon = " << theta_phonon
	     << " phi_phonon = " << phi_phonon << " q = " << q << G4endl;
    }
    
    // Sanity check for phonon production: should be forward, like Cherenkov
    if (theta_phonon>acos(kSound/kmag) || theta_phonon>halfpi) {
      if (verboseLevel) {
	G4cerr << GetProcessName() << " TRY AGAIN: theta_phonon "
	       << theta_phonon << " exceeds cone angle " << acos(kSound/kmag)
	       << G4endl;
      }

      continue;			// Try again
    }
    
    // Generate phonon momentum vector
    qvec = q*kdir;
    qvec.rotate(kdir.orthogonal(), theta_phonon);
    qvec.rotate(kdir, phi_phonon);
    
    if (verboseLevel > 1) {
      G4cout << " qvec = " << qvec << G4endl
	     << " ktrk.qvec = " << ktrk.dot(qvec)/(kmag*q)
	     << " ktr-qvec angle " << acos(ktrk.dot(qvec)/(kmag*q))
	     << G4endl;
    }
    
    // Get recoil wavevector (in HV frame), convert to new local momentum
    k_recoil = ktrk - qvec;
    if (IsHole()) {
      precoil = k_recoil * hbarc;
      Erecoil = precoil.mag2() / (2.*lat->GetHoleMass());
    } else {
      precoil = lat->MapK_HVtoP(iValley, k_recoil);
      Erecoil = lat->MapPtoEkin(iValley, precoil);
    }

    if (verboseLevel > 1) {
      G4cout << "k_recoil = " << k_recoil << " k_mag " << k_recoil.mag()
	     << G4endl;
      
      if (IsElectron()) {
	G4ThreeVector kvalley = theLattice->MapK_HVtoK_valley(iValley, k_recoil);
	G4ThreeVector pvalley = kvalley * hbarc;
	
	G4cout << "k_valley = " << kvalley << " k_mag " << kvalley.mag() << G4endl
	       << "p_valley = " << pvalley << " p_mag " << pvalley.mag() << G4endl;
      }
      
      G4cout << "p_recoil = " << precoil << " p_mag " << precoil.mag() << G4endl
	     << "E_recoil = " << Erecoil/eV << " eV" << G4endl;
    }
    
    // Sanity check for phonon production: can't exceed charge's energy
    Ephonon = MakePhononEnergy(qvec.mag());
    if (Ephonon >= GetKineticEnergy(aTrack)) {
      if (verboseLevel) {
	G4cerr << GetProcessName() << " TRY AGAIN: Ephonon "
	       << Ephonon/eV << " eV exceeds " << trkName << " energy "
	       << GetKineticEnergy(aTrack)/eV << " eV" << G4endl;
      }

      continue;			// Try again
    }

    // Sanity check for electrons: recoil energy must be smaller
    if (IsElectron()) {
      /***** Leave code in place for future tuning of IV scattering
      // See if momentum should be kicked into new valley
      newValley = FindNearestValley(precoil);

      if (verboseLevel && iValley != newValley) {
	G4double ErecNew = lat->MapPtoEkin(newValley, precoil);

	G4cout << " MOVED valley " << iValley << " to " << newValley
	       << " Erecoil : original " << Erecoil/eV << " new "
	       << ErecNew/eV << " eV" << G4endl;

	Erecoil = ErecNew;

	// Rescale phonon to enforce energy conservation; breaks momentum
	G4double Escale = (GetKineticEnergy(aTrack)-Erecoil) / Ephonon;
	qvec *= Escale;
	Ephonon *= Escale;
      }
      *****/

      G4double Efinal = Erecoil + Ephonon;
      G4double DeltaE = Efinal - GetKineticEnergy(aTrack);

      // If too much outgoing energy, throw again
      if (DeltaE > 1e-9) {
	if (verboseLevel) {
	  G4cerr << " TRY AGAIN: E(recoil+phonon) " << Efinal/eV << " eV"
		 << " too large by " << DeltaE/eV << G4endl;
	}
	
	continue;
      }

      // If too little outgoing energy, rescale momentum to balance
      if (DeltaE < -1e-9) {
	if (verboseLevel) {
	  G4cerr << " ECONS: E(recoil+phonon) " << Efinal/eV << " eV"
		 << " too small by " << DeltaE/eV << G4endl;
	}

	
      }
    }

    goodThrow = true;		// Nothing failed, get out of loop
  }	// while (goodThrow...)

  if (!goodThrow) {
    G4cerr << GetProcessName() << " ERROR: Unable to generate phonon after "
	   << iThrow << " attempts" << G4endl;
    return &aParticleChange;	// Unable to generate phonon
  } else if (verboseLevel && iThrow>1) {
    G4cout << GetProcessName() << " " << trkName << " phonon required "
	   << iThrow << " attempts" << G4endl;
  }

  // Report phonon emission results
  if (verboseLevel > 1) {
    G4cout << "q = " << q << " |qvec| = " << qvec.mag()
	   << " Ephonon = " << Ephonon
           << "\nk_recoil(HV) = " << k_recoil
           << " k_recoil(HV)-mag = " << k_recoil.mag()
	   << " newValley = " << newValley
           << G4endl;
  }

#ifdef G4CMP_DEBUG
  if (output.good()) {
    output << aTrack.GetTrackID() << "," << trkName << ","
	   << aTrack.GetWeight() << "," << GetKineticEnergy(aTrack)/eV << ","
	   << GetLocalMomentum(aTrack).mag()/eV << "," << kmag << ","
	   << theta_phonon << "," << Ephonon/eV << ",";
  }
#endif

  // Create real phonon to be propagated, with random polarization
  // If phonon is not created, register the energy as deposited
  G4double weight =
    G4CMP::ChoosePhononWeight(G4CMPConfigManager::GetLukeSampling());
  if (weight > 0.) {
    G4Track* phonon = G4CMP::CreatePhonon(aTrack,
					  G4PhononPolarization::UNKNOWN,
                                          qvec, Ephonon,
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

  RotateToGlobalDirection(precoil);	// Update track in world coordinates
  FillParticleChange(newValley, precoil);

#ifdef G4CMP_DEBUG
  if (output.good()) {
    output << aTrack.GetWeight()*weight << "," << k_recoil.mag() << ","
	   << Erecoil/eV << "," << precoil.mag()/eV << std::endl;
  }
#endif

  ClearNumberOfInteractionLengthLeft();
  return &aParticleChange;
}

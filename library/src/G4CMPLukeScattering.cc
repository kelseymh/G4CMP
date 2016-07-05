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

#include "G4CMPLukeScattering.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4Field.hh"
#include "G4FieldManager.hh"
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
#include <functional>


// Constructor and destructor

G4CMPLukeScattering::G4CMPLukeScattering()
  : G4CMPVDriftProcess("G4CMPLukeScattering", fLukeScattering) {
  LoadDataFromLUT(G4CMPConfigManager::GetLatticeDir() +
		  "/Ge/luke_scatter_times_elec.dat", ElecGPILData);
  LoadDataFromLUT(G4CMPConfigManager::GetLatticeDir() +
		  "/Ge/luke_scatter_times_hole.dat", HoleGPILData);

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
G4CMPLukeScattering::PostStepGetPhysicalInteractionLength(
                        const G4Track& aTrack,
                        G4double, G4ForceCondition* condition) {
  *condition = NotForced;

  G4double velocity = GetVelocity(aTrack);
  G4double mass = 0.;
  PDFDataTensor* GPILData = nullptr;
  G4ThreeVector k0;
  if (GetCurrentParticle() == G4CMPDriftElectron::Definition()) {
    mass = theLattice->GetElectronMass();
    GPILData = &ElecGPILData;
    k0 = theLattice->MapV_elToK_HV(GetValleyIndex(aTrack),
                                   GetGlobalVelocityVector(aTrack));
  } else if (GetCurrentParticle() == G4CMPDriftHole::Definition()) {
    mass = theLattice->GetHoleMass();
    GPILData = &HoleGPILData;
    k0 = GetLocalWaveVector(aTrack);
  } else {
    G4Exception("G4CMPLukeScattering::PostStepGetPhysicalInteractionLength",
                "G4CMPLuke003", EventMustBeAborted,
                "Invalid particle for LukeScatter process");
    return 0.;
  }
  //FIXME: This should be a global wave vector, not local?
  G4double kSound = CalculateKSound(mass);
  G4double mach = k0.mag()/kSound;

  if (verboseLevel > 1) {
    G4cout << "LukeScattering v = " << velocity/m*s << " kmag = " << k0.mag()*m
	   << G4endl;
  }

  G4FieldManager* fMan =
    aTrack.GetVolume()->GetLogicalVolume()->GetFieldManager();

  G4ThreeVector E0(0., 0., 0.);
  if (fMan && fMan->DoesFieldExist()) {
    const G4Field* field = fMan->GetDetectorField();
    //FIXME: This should be global position?
    G4ThreeVector x0 = GetLocalPosition(aTrack);
    G4double fieldValue[6];
    field->GetFieldValue(&x0[0], fieldValue);

    E0 = G4ThreeVector(fieldValue[3], fieldValue[4], fieldValue[5]);

    if (GetCurrentParticle() == G4CMPDriftElectron::Definition()) {
      E0 = theLattice->MapPtoK_HV(GetValleyIndex(aTrack), E0) * hbarc;
    } else if (GetCurrentParticle() == G4CMPDriftHole::Definition()) {
    } else {
      G4Exception("G4CMPLukeScattering::PostStepDoIt", "G4CMPLuke006",
                  EventMustBeAborted, "Invalid particle for LukeScatter process");
    }
  }

  G4double theta = 0;
  if (k0.mag() && E0.mag()) {
    G4double arg = k0*E0/(k0.mag()*E0.mag());
    if (arg > 1.) arg = 1.;
    if (arg < -1.) arg = -1.;
    theta = acos(arg);
  }

  size_t thetabin = 0;
  if (theta > THETAMAX) {
    thetabin = THETASIZE - 1;
  } else if (THETASIZE > 1 && THETAMAX-THETAMIN > 0) {
    thetabin = round(theta-THETAMIN)/((THETAMAX-THETAMIN)/(THETASIZE-1));
  }

  size_t machbin = 0;
  if (mach > MACHMAX) {
    machbin = MACHSIZE - 1;
  } else if (MACHSIZE > 1 && MACHMAX-MACHMIN > 0) {
    machbin = round(mach-MACHMIN)/((MACHMAX-MACHMIN)/(MACHSIZE-1));
  }

  size_t ebin = 0;
  if (E0.mag() > EMAX) {
    ebin = ESIZE - 1;
  } else if (ESIZE > 1 && EMAX-EMIN > 0) {
    ebin = round(E0.mag()-EMIN)/((EMAX-EMIN)/(ESIZE-1));
  }

  if (thetabin >= THETASIZE || ebin >= ESIZE || machbin >= MACHSIZE) {
    G4Exception("G4CMPLukeScattering::PostStepGetPhysicalInteractionLength",
                "G4CMPLuke007", EventMustBeAborted,
                "Invalid index for Luke phonon scatter look-up table.");
  }

  G4double code = (*GPILData)[ebin][thetabin][machbin][0];
  G4double mean = (*GPILData)[ebin][thetabin][machbin][1];
  G4double std  = (*GPILData)[ebin][thetabin][machbin][2];

  std::function<G4double(G4double, G4double)> shootRand;
  if (code == 0) {
    shootRand = [](G4double avg, G4double stddev) {
                  return G4RandGauss::shoot(avg, stddev); };
  } else if (code == 1) {
    shootRand = [](G4double avg, G4double /*stddev*/) {
                  return G4RandExponential::shoot(avg); };
  } else if (code == 2) {
    if (verboseLevel > 1) {
      G4cout << "LukeScattering PIL = DBL_MAX" << G4endl;
    }
    return DBL_MAX;
  } else {
    G4ExceptionDescription
        desc("LukeScattering lookup table invalid code for GPIL");
    G4Exception("G4CMPLukeScattering::GetPhysicalInteractionLength",
                "G4CMPLuke002", EventMustBeAborted, desc);
  }

  G4double pil = 0;
  G4double dt = 0;
  G4ThreeVector dx(0.);
  G4int count = 0;
  while (pil <= 0.) {
    dt = shootRand(mean, std);
    dx = (hbar_Planck * k0 * dt +
           0.5 * aTrack.GetDynamicParticle()->GetCharge() * E0 * dt * dt)
          / mass;

    if (GetCurrentParticle() == G4CMPDriftElectron::Definition()) {
      dx *= theLattice->GetSqrtInvTensor();
      dx.transform(GetValley(aTrack).inverse());
    } else if (GetCurrentParticle() == G4CMPDriftHole::Definition()) {
    } else {
      G4Exception("G4CMPLukeScattering::PostStepDoIt", "G4CMPLuke006",
                  EventMustBeAborted, "Invalid particle for LukeScatter process");
    }

    pil = dx.mag();

    if (++count >= 1000) {
      G4ExceptionDescription
          desc("Spending too long trying to generate a PIL");
      G4Exception("G4CMPLukeScattering::GetPhysicalInteractionLength",
                  "G4CMPLuke004", EventMustBeAborted, desc);
    }
  }

  if (verboseLevel > 1) {
    G4cout << "LukeScattering PIL = " << pil << G4endl;
  }

  return pil;
}


G4VParticleChange* G4CMPLukeScattering::PostStepDoIt(const G4Track& aTrack,
                                                     const G4Step& aStep) {
  aParticleChange.Initialize(aTrack); 
  G4StepPoint* postStepPoint = aStep.GetPostStepPoint();

  G4double mass = 0.;
  if (GetCurrentParticle() == G4CMPDriftElectron::Definition()) {
    mass = theLattice->GetElectronMass();
  } else if (GetCurrentParticle() == G4CMPDriftHole::Definition()) {
    mass = theLattice->GetHoleMass();
  } else {
    G4Exception("G4CMPLukeScattering::PostStepDoIt", "G4CMPLuke005",
                EventMustBeAborted, "Invalid particle for LukeScatter process");
    return &aParticleChange;
  }
  
  G4ThreeVector ktrk = GetLocalWaveVector(aTrack);
  if (GetCurrentParticle() == G4CMPDriftElectron::Definition()) {
    ktrk = theLattice->MapV_elToK_HV(GetValleyIndex(aTrack),
                                     GetGlobalVelocityVector(aTrack));
  }
  G4double kmag = ktrk.mag();
  G4double kSound = CalculateKSound(mass);

  // GPIL returns a length that is too short for very few cases because
  // the PDF in the "Gaussian" case is not truly Gaussian at the low tail.
  if (kmag <= kSound) {
    return &aParticleChange;
  }

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
  if (genLuke > 0. && G4UniformRand() < genLuke) {
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

void G4CMPLukeScattering::LoadDataFromLUT(const G4String &filename,
					  PDFDataTensor& GPILData) {
  std::ifstream LUT(filename);

  std::string temp;
  std::getline(LUT,temp); // First header

  if (LUT.good()) {
    LUT >> ESIZE >> EMIN >> EMAX
        >> THETASIZE >> THETAMIN >> THETAMAX
        >> MACHSIZE >> MACHMIN >> MACHMAX;
        EMIN *= volt/m;
        EMAX *= volt/m;
        THETAMIN *= rad;
        THETAMAX *= rad;
  } else {
    G4ExceptionDescription desc("Couldn't find and/or open luke_mfp.dat in "
                                + G4CMPConfigManager::GetLatticeDir());
    G4Exception("G4CMPLukeScattering::LoadDataFromLUT()",
                "G4CMPLuke001", FatalException, desc);
  }

  std::getline(LUT,temp); // Newline char
  std::getline(LUT,temp); // Second header

  G4double code, mean, std;
  G4String fitType;
  GPILData.resize(ESIZE,
		  PDFDataMatrix(THETASIZE,
				PDFDataRow(MACHSIZE,
					   PDFParamTuple())));
  for (size_t i=0; i<ESIZE; ++i) {
    for (size_t j=0; j<THETASIZE; ++j) {
      for (size_t k=0; k<MACHSIZE; ++k) {
        LUT >> fitType >> mean >> std;

        if (fitType == "gaus") code = 0;
        else if (fitType == "exp") code = 1;
        else if (fitType == "flat") code = 2;
        else G4cout << "BAD" << G4endl;

        GPILData[i][j][k][0] = code;
        GPILData[i][j][k][1] = mean * ns;
        GPILData[i][j][k][2] = std * ns;
      }
    }
  }
}

G4double G4CMPLukeScattering::CalculateKSound(G4double mass) {
  G4double vSound = theLattice->GetSoundSpeed();
  return vSound*mass/hbar_Planck;
}

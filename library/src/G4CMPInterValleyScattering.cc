/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
//
// 20140324  Drop hard-coded IV scattering parameters; get from lattice
// 20140324  Restore Z-axis mass tensor
// 20140331  Add required process subtype code
// 20140418  Drop local valley transforms, use lattice functions instead
// 20140429  Recompute kinematics relative to new valley
// 20140908  Allow IV scatter to change momentum by conserving energy
// 20150109  Revert IV scattering to preserve momentum
// 20150112  Follow renaming of "SetNewKinematics" to FillParticleChange
// 20150122  Use verboseLevel instead of compiler flag for debugging
// 20160601  Must apply lattice rotation before valley.
// 20161004  Change valley selection function to avoid null choice
// 20161114  Use G4CMPDriftTrackInfo
// 20170602  Use G4CMPUtils for particle identity checks

#include "G4CMPInterValleyScattering.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftTrackInfo.hh"
#include "G4CMPFieldManager.hh"
#include "G4CMPGeometryUtils.hh"
#include "G4CMPTrackUtils.hh"
#include "G4CMPUtils.hh"
#include "G4Field.hh"
#include "G4FieldManager.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalVolume.hh"
#include "G4PhononPolarization.hh"
#include "G4PhysicalConstants.hh"
#include "G4RandomDirection.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4CMPUtils.hh"
#include "G4VParticleChange.hh"
#include "G4VPhysicalVolume.hh"
#include "Randomize.hh"
#include "math.h"

G4CMPInterValleyScattering::G4CMPInterValleyScattering()
  : G4CMPVDriftProcess("G4CMPInterValleyScattering", fInterValleyScattering) {;}

G4CMPInterValleyScattering::~G4CMPInterValleyScattering() {;}


<<<<<<< HEAD
G4bool 
G4CMPInterValleyScattering::IsApplicable(const G4ParticleDefinition& aPD) {
  return G4CMP::IsElectron(&aPD);
}


G4double 
G4CMPInterValleyScattering::GetMeanFreePath(const G4Track& aTrack,
					    G4double,
					    G4ForceCondition* condition) {
  *condition = NotForced;
=======
G4double
G4CMPInterValleyScattering::GetMeanFreePath(const G4Track& aTrack, G4double,
                                     G4ForceCondition* condition) {
  *condition = Forced;		// In order to recompute MFP after TimeStepper

  G4double kmag = 0.;
  if (IsElectron(&aTrack)) {
    kmag = theLattice->MapV_elToK_HV(GetValleyIndex(aTrack),
                                     GetLocalVelocityVector(aTrack)).mag();
  }

  G4CMPTrackInformation* trackInfo = GetTrackInfo(aTrack);
  G4double kSound = CalculateKSound(trackInfo);
  G4double velocity = GetVelocity(aTrack);
  G4double Eelectron = GetKineticEnergy(aTrack);
>>>>>>> 58a775f003c357a55b7d156c3681999ac1de2ff3

  G4FieldManager* fMan =
    aTrack.GetVolume()->GetLogicalVolume()->GetFieldManager();
  //If there is no field, there is no IV scattering... but then there
  //is no e-h transport either...
  if (!fMan || !fMan->DoesFieldExist()) return DBL_MAX;
  if (kmag <= kSound) return DBL_MAX;

  // these numbers from V. Aubry-Fortuna PACS 72.10.Di, 72.20.Fr
  G4double D1 = 3 * pow(10,8)/eV*cm; // deformation potential for first transition
  G4double E1 = 0.0276/eV; // phonon energy for first transition
  G4double D2 = 2 * pow(10, 7)/eV*cm; // deformation potential for second transition
  G4double E2 = 0.0103/eV; // phonon energy for second transition

  //FIXME: impurity scattering place-holder
  G4double impurityRate = 0;

  G4double rate = GetIVSRate(D1, E1, Eelectron) + GetIVSRate(D2, E2, Eelectron)
                  + impurityRate;
  G4double l0 = rate/velocity;

  // Time step corresponding to Mach number (avg. time between radiations)
  G4double dtau = ChargeCarrierTimeStep(kmag/kSound, l0);
  G4double mfp = dtau * velocity;

  return mfp;
}

G4VParticleChange* G4CMPInterValleyScattering::PostStepDoIt(const G4Track& aTrack,
                                                            const G4Step& aStep) {
  aParticleChange.Initialize(aTrack);
  G4StepPoint* postStepPoint = aStep.GetPostStepPoint();

  if (postStepPoint->GetStepStatus()==fGeomBoundary) {
      return &aParticleChange;
  }

  G4ThreeVector ktrk(0.);
  if (IsElectron(&aTrack)) {
    ktrk = theLattice->MapV_elToK_HV(GetValleyIndex(aTrack),
                                     GetLocalVelocityVector(aTrack));
  }

  G4double kmag = ktrk.mag();
  G4double kSound = CalculateKSound(GetTrackInfo(aTrack));

  // Sanity check: this should have been done in MFP already
  if (kmag <= kSound) return &aParticleChange;

<<<<<<< HEAD
  // picking a new valley at random if IV-scattering process was triggered
  valley = ChangeValley(valley);
  G4CMP::GetTrackInfo<G4CMPDriftTrackInfo>(aTrack)->SetValleyIndex(valley);
=======
  G4double theta_phonon = MakePhononTheta(kmag, kSound);
  G4double phi_phonon   = G4UniformRand()*twopi;
  G4double q = 2*(kmag*cos(theta_phonon)-kSound);
>>>>>>> 58a775f003c357a55b7d156c3681999ac1de2ff3

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

  G4double Ephonon;
  G4double num = G4UniformRand();
  if (num < 1/3)
    Ephonon = 0.0276/eV;
  else if (num > 1/3)
    Ephonon = 0.0103/eV;

  // Get recoil wavevector, convert to new momentum
  G4ThreeVector k_recoil = ktrk - qvec;

  // Create real phonon to be propagated, with random polarization
  // If phonon is not created, register the energy as deposited
  G4double weight = G4CMP::ChoosePhononWeight();
  if (weight > 0.) {
    MakeGlobalPhononK(qvec);  		// Convert phonon vector to real space

    G4Track* phonon = CreatePhonon(G4PhononPolarization::UNKNOWN,qvec,Ephonon);
    phonon->SetWeight(weight);

    aParticleChange.SetNumberOfSecondaries(1);
    aParticleChange.AddSecondary(phonon);

    //deciding which valley to scatter into
    G4int valley = GetValleyIndex(aTrack);
    G4int newValley = 0;
    if (Ephonon == 0.0103/eV) {
      G4double whichValley = G4UniformRand();
      if (whichValley < 0.5)
        newValley = (valley + 1) % 4;
      else if (whichValley > 0.5)
        newValley = (valley + 3) % 4;
      GetTrackInfo(aTrack)->SetValleyIndex(newValley);
    }
    if (Ephonon == 0.0276/eV) {
        newValley = (valley + 2) % 4;
        GetTrackInfo(aTrack)->SetValleyIndex(newValley);
    }

    MakeGlobalRecoil(k_recoil);
    FillParticleChange(valley, k_recoil);
  } else {
    aParticleChange.ProposeNonIonizingEnergyDeposit(Ephonon);

    MakeGlobalRecoil(k_recoil);
    FillParticleChange(GetValleyIndex(aTrack), k_recoil);
  }

  ResetNumberOfInteractionLengthLeft();
  return &aParticleChange;
}
<<<<<<< HEAD
=======


G4double G4CMPInterValleyScattering::GetIVSRate(const G4double defpot,
                                                      const G4double Ephonon,
                                                      const G4double Eelectron) {
  // calculates an intervalley scattering rate based on
  // V. Aubry-Fortuna et. al.
  // Solid-State Electronics 49 (2005) 1320-1329
  G4int Ziv = 3; // number of possible final valleys
  G4double mLong = theLattice->GetMassTensor()[1][1];
  G4double mTrans = theLattice->GetMassTensor()[2][2];
  G4double mDOS = cbrt(mTrans*mTrans*mLong);
  G4double density = theLattice->GetLattice()->GetDensity();
  G4double DeltaEiv = 0; // energy difference between valleys
  G4double alpha = 0.3*eV; // non-parabolicity as per V. Aubry-Fortuna PACS 72.10.Di, 72.20.Fr

  G4double IVSrate = Ziv/sqrt(2)/pi * pow(mDOS, 1.5)*defpot*defpot/(hbar_Planck*hbar_Planck)
                     /density/Ephonon * sqrt(Eelectron + Ephonon + DeltaEiv)
                     * sqrt(1+alpha*(Eelectron + Ephonon + DeltaEiv))
                     * (1+2*alpha*(Eelectron + Ephonon + DeltaEiv));

  return IVSrate;
}

G4double
G4CMPInterValleyScattering::CalculateKSound(const G4CMPTrackInformation* trackInfo) {
  return theLattice->GetSoundSpeed()*trackInfo->GetEffectiveMass()/hbar_Planck;
}

G4bool G4CMPInterValleyScattering::IsApplicable(const G4ParticleDefinition& aPD)
{
  return (&aPD==G4CMPDriftElectron::Definition());
}
>>>>>>> 58a775f003c357a55b7d156c3681999ac1de2ff3

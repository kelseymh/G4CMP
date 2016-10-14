/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

//
/// \file library/src/G4CMPProcessUtils.cc
/// \brief Implementation of the G4CMPProcessUtils class
///   Provides useful general functions to fetch and store lattice, access
///   and apply lattice parameters for phonons and charge carriers.
///
///   Use via multiple inheritance with concrete or base process classes
//
// $Id$
//
// 20140321  Move lattice-based placement transformations here, via Touchable
// 20140407  Add functions for phonon generation in Luke scattering
// 20140412  Add manual configuration options
// 20140509  Add ChoosePolarization() which uses DOS values from lattice
// 20141216  Set velocity "by hand" for secondary electrons
// 20150109  Use G4CMP_SET_ELECTRON_MASS to enable dynamic mass, velocity set
// 20150112  Add GetCurrentValley() function to get valley of current track,
//	     allow GetValley functions to treat holes, returning -1
// 20150309  Add Create*() functions which take position and energy arguments
//	     (for use with AlongStepDoIt() actions).
// 20150310  Fix CreateChargeCarrier to use momentum unit vector
// 20160610  Return regular (NOT Herring-Vogt) wave vector for electrons
// 20160809  BUG FIX:  th_phonon==0 is fine for computing energy.
// 20160825  Add assignment operators for cross-process configuration;
//	     move track identification functions to G4CMPUtils.
// 20160829  Drop G4CMP_SET_ELECTRON_MASS code blocks; not physical
// 20160906  Make GetSurfaceNormal() const.
// 20161004  Add new ChangeValley() function to avoid null selection

#include "G4CMPProcessUtils.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4CMPTrackInformation.hh"
#include "G4CMPUtils.hh"
#include "G4AffineTransform.hh"
#include "G4DynamicParticle.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParallelWorldProcess.hh"
#include "G4PhononLong.hh"
#include "G4PhononPolarization.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononTransSlow.hh"
#include "G4PhysicalConstants.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4TransportationManager.hh"
#include "G4VTouchable.hh"
#include "Randomize.hh"

#include "G4GeometryTolerance.hh"
#include "G4VSolid.hh"

// Constructor and destructor

G4CMPProcessUtils::G4CMPProcessUtils()
  : theLattice(nullptr), fPhysicsModelID(0), currentTrack(nullptr),
    currentVolume(nullptr) {
  fPhysicsModelID = G4PhysicsModelCatalog::Register("G4CMP process");
}

G4CMPProcessUtils::~G4CMPProcessUtils() {;}


// Assignment operators allow dependent configuration

G4CMPProcessUtils::G4CMPProcessUtils(G4CMPProcessUtils& right)
  : theLattice(right.theLattice), fPhysicsModelID(right.fPhysicsModelID),
    currentTrack(right.currentTrack), currentVolume(right.currentVolume),
    fLocalToGlobal(right.fLocalToGlobal),
    fGlobalToLocal(right.fGlobalToLocal) {;}

G4CMPProcessUtils&
G4CMPProcessUtils::operator=(const G4CMPProcessUtils& right) {
  if (this != &right) {			// Avoid unnecessary work
    theLattice      = right.theLattice;
    fPhysicsModelID = right.fPhysicsModelID;
    currentTrack    = right.currentTrack;
    currentVolume   = right.currentVolume;
    fLocalToGlobal  = right.fLocalToGlobal;
    fGlobalToLocal  = right.fGlobalToLocal;
  }
  return *this;
}


// Initialization for current track

void G4CMPProcessUtils::LoadDataForTrack(const G4Track* track) {
  currentTrack = track;
  currentVolume = track->GetVolume();
//  G4cout << "New track:" << G4endl;
//  G4cout << track->GetDefinition()->GetParticleName() << G4endl;
//  G4cout << track->GetPosition() << G4endl;
//  G4cout << track->GetMomentumDirection() << G4endl;
//  if (aTrack->GetCreatorProcess()) {
//  G4cout << "Creator Process: " << aTrack->GetCreatorProcess()->GetProcessName() << G4endl;
//  }
//  if (aTrack->GetMaterial()) {
//  G4cout << "Material: " << aTrack->GetMaterial()->GetName() << G4endl;
//  }
//  if (track->GetVolume()) {
//  if (track->GetCreatorProcess())
//  G4cout << "Parent Proc: " << track->GetCreatorProcess()->GetProcessName() << G4endl;
//  G4cout << "Log vol at vert: " << track->GetLogicalVolumeAtVertex()->GetName() << G4endl;
//  G4cout << "Origin touchable: " << track->GetOriginTouchable()->GetVolume(0)->GetName() << G4endl;
//  G4cout << "Touchable: " << track->GetTouchable()->GetVolume(0)->GetName() << G4endl;
//  G4cout << "NextTouchable: " << track->GetNextTouchable()->GetVolume(0)->GetName() << G4endl;
//  G4cout << "Volume: " << track->GetVolume()->GetName() << G4endl;
//  }

  // WARNING!  This assumes track starts and ends in one single volume!
  SetTransforms(track->GetTouchable());
  FindLattice(track->GetVolume());

  if (!theLattice) {
    G4Exception("G4CMPProcessUtils::LoadDataForTrack", "Utils001",
		EventMustBeAborted, ("No lattice found for volume"
				     + track->GetVolume()->GetName()).c_str());
    return;	// No lattice, no special actions possible
  }

  // Fill auxiliary track info not already recorded (needs non-const object)
  G4CMPTrackInformation* trackInfo = AttachTrackInfo(track);

  if (IsPhonon(track)) {
    // Set momentum direction using already provided wavevector
    G4ThreeVector kdir = trackInfo->GetPhononK();

    const G4ParticleDefinition* pd = track->GetParticleDefinition();
    G4Track* tmp_track = const_cast<G4Track*>(track);
    tmp_track->SetMomentumDirection(
      theLattice->MapKtoVDir(G4PhononPolarization::Get(pd), kdir));
  }

  if (IsElectron(track)) {
    trackInfo->SetScatterLength(theLattice->GetElectronScatter());
    trackInfo->SetEffectiveMass(theLattice->GetElectronMass());
    if (trackInfo->GetValleyIndex() < 0)
      trackInfo->SetValleyIndex(ChooseValley());
  }

  if (IsHole(track)) {
    trackInfo->SetScatterLength(theLattice->GetHoleScatter());
    trackInfo->SetEffectiveMass(theLattice->GetHoleMass());
    trackInfo->SetValleyIndex(-1);		// Holes don't have valleys
  }
}

// Create new info object or return existing one for track

G4CMPTrackInformation*
G4CMPProcessUtils::AttachTrackInfo(const G4Track* track) const {
  if (track == 0) return nullptr;		// Must have valid track

  G4CMPTrackInformation* trkInfo = GetTrackInfo(track);
  if (!trkInfo) {
    trkInfo = new G4CMPTrackInformation;
    track->SetAuxiliaryTrackInformation(fPhysicsModelID, trkInfo);
  }

  return trkInfo;
}


// Identify track type to simplify some conditionals

G4bool G4CMPProcessUtils::IsPhonon(const G4Track* track) const {
  return G4CMP::IsPhonon(track);
}

G4bool G4CMPProcessUtils::IsElectron(const G4Track* track) const {
  return G4CMP::IsElectron(track);
}

G4bool G4CMPProcessUtils::IsHole(const G4Track* track) const {
  return G4CMP::IsHole(track);
}

G4bool G4CMPProcessUtils::IsChargeCarrier(const G4Track* track) const {
  return G4CMP::IsChargeCarrier(track);
}



// Fetch lattice for current track, use in subsequent steps

void G4CMPProcessUtils::FindLattice(const G4VPhysicalVolume* volume) {
  currentVolume = volume;		// 

  G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();
  theLattice = LM->GetLattice(volume);

  if (!theLattice) {
    G4cerr << "WARNING: No lattice for volume " << volume->GetName() << G4endl;
  }
}

// Configure orientation matrices for current track

void G4CMPProcessUtils::SetTransforms(const G4VTouchable* touchable) {
  if (!touchable) {			// Null pointer defaults to identity
    fLocalToGlobal = fGlobalToLocal = G4AffineTransform();
    return;
  }

  SetTransforms(touchable->GetRotation(), touchable->GetTranslation());
}

void G4CMPProcessUtils::SetTransforms(const G4RotationMatrix* rot,
				      const G4ThreeVector& trans) {
  fLocalToGlobal = G4AffineTransform(rot, trans);
  fGlobalToLocal = fLocalToGlobal.Inverse();
}

// Delete current configuration before new track starts

void G4CMPProcessUtils::ReleaseTrack() {
  SetTransforms(nullptr);
  currentTrack = nullptr;
  theLattice = nullptr;
}

G4ThreeVector G4CMPProcessUtils::GetSurfaceNormal(const G4Step& aStep) const {
  // Get outward normal using G4Navigator method (more reliable than G4VSolid)
  G4int navID = G4ParallelWorldProcess::GetHypNavigatorID();
  std::vector<G4Navigator*>::iterator iNav =
    G4TransportationManager::GetTransportationManager()->GetActiveNavigatorsIterator();

  G4bool goodNorm;
  G4ThreeVector surfNorm = iNav[navID]->GetGlobalExitNormal(
                                      aStep.GetPostStepPoint()->GetPosition(),
                                      &goodNorm);

  // FIXME:  Sometimes G4Navigator fails, but still returns "good"
  if (!goodNorm || surfNorm.mag()<0.99) {
    G4VPhysicalVolume* thePrePV = aStep.GetPreStepPoint()->GetPhysicalVolume();
    G4VPhysicalVolume* thePostPV = aStep.GetPostStepPoint()->GetPhysicalVolume();
    G4Exception("G4CMPProcessUtils::GetSurfaceNormal", "Boundary001",
                EventMustBeAborted, ("Can't get normal vector of surface between " +
                                    thePrePV->GetName() + " and " +
                                    thePostPV->GetName()+ ".").c_str());
  }
  return surfNorm;
}

// Access track position and momentum in local coordinates
G4ThreeVector G4CMPProcessUtils::GetLocalPosition(const G4Track& track) const {
  return GetLocalPosition(track.GetPosition());
}

void G4CMPProcessUtils::GetLocalPosition(const G4Track& track,
					 G4double pos[3]) const {
  G4ThreeVector tpos = GetLocalPosition(track);
  pos[0] = tpos.x();
  pos[1] = tpos.y();
  pos[2] = tpos.z();
}

G4ThreeVector G4CMPProcessUtils::GetLocalMomentum(const G4Track& track) const {
  if (IsElectron(&track)) {
    return theLattice->MapV_elToP(GetValleyIndex(track),
                                  GetLocalVelocityVector(track));
  } else if (IsHole(&track)) {
    return GetLocalDirection(track.GetMomentum());
  } else {
    G4Exception("G4CMPProcessUtils::GetLocalMomentum()", "DriftProcess001",
                EventMustBeAborted, "Unknown charge carrier");
    return G4ThreeVector();
  }
}

void G4CMPProcessUtils::GetLocalMomentum(const G4Track& track, 
					 G4double mom[3]) const {
  G4ThreeVector tmom = GetLocalMomentum(track);
  mom[0] = tmom.x();
  mom[1] = tmom.y();
  mom[2] = tmom.z();
}

G4ThreeVector G4CMPProcessUtils::GetLocalVelocityVector(const G4Track& track) const {
  G4ThreeVector vel = track.CalculateVelocity() * track.GetMomentumDirection();
  return GetLocalDirection(vel);
}

void G4CMPProcessUtils::GetLocalVelocityVector(const G4Track &track,
                                               G4double vel[]) const {
  G4ThreeVector v_local = GetLocalVelocityVector(track);
  vel[0] = v_local.x();
  vel[1] = v_local.y();
  vel[2] = v_local.z();
}

G4ThreeVector G4CMPProcessUtils::GetLocalWaveVector(const G4Track& track) const {
  if (IsChargeCarrier(&track)) {
    return GetLocalMomentum(track) / hbarc;
  } else if (IsPhonon(&track)) {
    return GetTrackInfo(track)->GetPhononK();
  } else {
    G4Exception("G4CMPProcessUtils::GetLocalWaveVector", "DriftProcess002",
                EventMustBeAborted, "Unknown charge carrier");
    return G4ThreeVector();
  }
}

// Access track position and momentum in global coordinates
G4ThreeVector G4CMPProcessUtils::GetGlobalPosition(const G4Track& track) const {
  return track.GetPosition();
}

void G4CMPProcessUtils::GetGlobalPosition(const G4Track& track,
           G4double pos[3]) const {
  G4ThreeVector tpos = GetGlobalPosition(track);
  pos[0] = tpos.x();
  pos[1] = tpos.y();
  pos[2] = tpos.z();
}

G4ThreeVector G4CMPProcessUtils::GetGlobalMomentum(const G4Track& track) const {
  if (IsElectron(&track)) {
    G4ThreeVector p = theLattice->MapV_elToP(GetValleyIndex(track),
                                             GetLocalVelocityVector(track));
    return GetGlobalDirection(p);
  } else if (IsHole(&track)) {
    return track.GetMomentum();
  } else {
    G4Exception("G4CMPProcessUtils::GetGlobalMomentum", "DriftProcess003",
                EventMustBeAborted, "Unknown charge carrier");
    return G4ThreeVector();
  }
}

void G4CMPProcessUtils::GetGlobalMomentum(const G4Track& track,
					  G4double mom[3]) const {
  G4ThreeVector tmom = GetGlobalMomentum(track);
  mom[0] = tmom.x();
  mom[1] = tmom.y();
  mom[2] = tmom.z();
}

G4ThreeVector G4CMPProcessUtils::GetGlobalVelocityVector(const G4Track& track) const {
  return track.CalculateVelocity() * track.GetMomentumDirection();
}

void G4CMPProcessUtils::GetGlobalVelocityVector(const G4Track &track, G4double vel[]) const {
  G4ThreeVector v_local = GetGlobalVelocityVector(track);
  vel[0] = v_local.x();
  vel[1] = v_local.y();
  vel[2] = v_local.z();
}

G4double G4CMPProcessUtils::GetKineticEnergy(const G4Track &track) const {
  if (IsElectron(&track)) {
    return theLattice->MapV_elToEkin(GetValleyIndex(track),
                                     GetLocalVelocityVector(track));
  } else if (IsHole(&track)) {
    return track.GetKineticEnergy();
  } else if (IsPhonon(&track)) {
    return track.GetKineticEnergy();
  } else {
    G4Exception("G4CMPProcessUtils::GetKineticEnergy", "G4CMPProcess004",
                EventMustBeAborted, "Unknown condensed matter particle");
    return 0.0;
  }
}

// Return particle type for currently active track [set in LoadDataForTrack()]

const G4ParticleDefinition* G4CMPProcessUtils::GetCurrentParticle() const {
  return (currentTrack ? currentTrack->GetParticleDefinition() : 0);
}


// Return auxiliary information for track (phonon, charge kinematics)

G4CMPTrackInformation* 
G4CMPProcessUtils::GetTrackInfo(const G4Track* track) const {
  if (!track) track = GetCurrentTrack();	// No argument, use current
  if (!track) return 0;

  return dynamic_cast<G4CMPTrackInformation*>
    (track->GetAuxiliaryTrackInformation(fPhysicsModelID));
}


// Access phonon particle-type/polarization indices

G4int G4CMPProcessUtils::GetPolarization(const G4Track& track) const {
  return G4PhononPolarization::Get(track.GetParticleDefinition());
}


// Generate random polarization from density of states

G4int G4CMPProcessUtils::ChoosePhononPolarization() const {
  return G4CMP::ChoosePhononPolarization(theLattice->GetLDOS(),
                                         theLattice->GetSTDOS(),
                                         theLattice->GetFTDOS());
}

void G4CMPProcessUtils::MakeLocalPhononK(G4ThreeVector& kphonon) const {
  if (IsElectron(GetCurrentTrack())) {
    kphonon = theLattice->MapK_HVtoK(GetValleyIndex(GetCurrentTrack()), kphonon);
  } else if (!IsHole(GetCurrentTrack())) {
    G4Exception("G4CMPProcessUtils::MakeGlobalPhonon", "DriftProcess005",
                EventMustBeAborted, "Unknown charge carrier");
  }
}

void G4CMPProcessUtils::MakeGlobalPhononK(G4ThreeVector& kphonon) const {
  MakeLocalPhononK(kphonon);
  RotateToGlobalDirection(kphonon);
}

void G4CMPProcessUtils::MakeGlobalRecoil(G4ThreeVector& kphonon) const {
  if (IsElectron(GetCurrentTrack())) {
    kphonon = theLattice->MapK_HVtoP(GetValleyIndex(GetCurrentTrack()),kphonon);
  } else if (IsHole(GetCurrentTrack())) {
    kphonon *= hbarc;
  } else {
    G4Exception("G4CMPProcessUtils::MakeGlobalPhonon", "DriftProcess006",
                EventMustBeAborted, "Unknown charge carrier");
  }
  RotateToGlobalDirection(kphonon);
}

// Compute a Lambertian distribution for reflected phonons

G4ThreeVector
G4CMPProcessUtils::LambertReflection(const G4ThreeVector& surfNorm) {
  G4double phi = 2.0*pi*G4UniformRand();
  G4double theta = acos(2.0*G4UniformRand() - 1.0) / 2.0;

  G4ThreeVector refl = -surfNorm;
  refl = refl.rotate(surfNorm.orthogonal(), theta);
  refl = refl.rotate(surfNorm, phi);
  return refl;
}

// Model Kaplan phonon-quasiparticle interactions in superconductor sensors

G4double G4CMPProcessUtils::KaplanPhononQP(G4double energy,
                                     G4MaterialPropertiesTable* prop,
                                     std::vector<G4double>& reflectedEnergies) {
  if (reflectedEnergies.size()>0)
    G4Exception("G4CMPProcessUtils::KaplanPhononQP()", "ProcessUtils007",
                JustWarning, "Passed a nonempty vector.");
  // Check that the MaterialPropertiesTable has everything we need. If it came
  // from a G4CMPSurfaceProperty, then it will be fine.
  if (!(prop->ConstPropertyExists("gapEnergy") &&
        prop->ConstPropertyExists("lowQPLimit") &&
        prop->ConstPropertyExists("phononLifetime") &&
        prop->ConstPropertyExists("phononLifetimeSlope") &&
        prop->ConstPropertyExists("vSound") &&
        prop->ConstPropertyExists("filmThickness"))) {
    G4Exception("G4CMPProcessUtils::KaplanPhononQP()", "ProcessUtils001",
                RunMustBeAborted,
                "Insufficient info in MaterialPropertiesTable.");
  }
  G4double gapEnergy     = prop->GetConstProperty("gapEnergy");
  G4double lowQPLimit    = prop->GetConstProperty("lowQPLimit");

  //For the phonon to not break a Cooper pair, it must go 2*thickness
  G4double frac = 2.0;

  G4double phononEscapeProb = CalcEscapeProbability(energy, frac, prop);

  G4double EDep = 0.;
  if (energy > 2.0*gapEnergy && G4UniformRand() > phononEscapeProb) {
    std::vector<G4double> qpEnergies;
    std::vector<G4double> phonEnergies{energy};
    while (qpEnergies.size()>0 || phonEnergies.size()>0) {
      if (phonEnergies.size()>0) {
        // Partition the phonons' energies into quasi-particles according to
        // a PDF defined in CalcQPEnergies().
        // NOTE: Both energy vectors mutate.
        EDep += CalcQPEnergies(gapEnergy, lowQPLimit, phonEnergies, qpEnergies);
      }
      if (qpEnergies.size()>0) {
        // Quasiparticles can also excite phonons.
        // NOTE: Both energy vectors mutate.
        EDep += CalcPhononEnergies(gapEnergy, lowQPLimit, phonEnergies, qpEnergies);
      }
      if (phonEnergies.size()>0) {
        // Some phonons will escape back into the crystal.
        // NOTE: Both energy vectors mutate.
        CalcReflectedPhononEnergies(prop, phonEnergies, reflectedEnergies);
      }
    }
  } else {
    reflectedEnergies.push_back(energy);
  }

  return EDep;
}

// Compute the probability of a phonon reentering the crystal

G4double G4CMPProcessUtils::CalcEscapeProbability(G4double energy,
                                              G4double thicknessFrac,
                                              G4MaterialPropertiesTable* prop) {
  G4double gapEnergy = prop->GetConstProperty("gapEnergy");
  G4double phononLifetime = prop->GetConstProperty("phononLifetime");
  G4double phononLifetimeSlope = prop->GetConstProperty("phononLifetimeSlope");
  G4double vSound = prop->GetConstProperty("vSound");
  G4double thickness = prop->GetConstProperty("filmThickness");

  G4double mfp = vSound * phononLifetime / (1. + phononLifetimeSlope * (
    energy/gapEnergy - 2.));
  return exp(-2.* thicknessFrac * thickness/mfp);
}

// Model the phonons breaking Cooper pairs into quasiparticles

G4double G4CMPProcessUtils::CalcQPEnergies(G4double gapEnergy,
                                           G4double lowQPLimit,
                                           std::vector<G4double>& phonEnergies,
                                           std::vector<G4double>& qpEnergies) {
  // Each phonon gives all of its energy to the qp pair it breaks.
  G4double EDep = 0.;
  for (G4double E: phonEnergies) {
    G4double qpE = QPEnergyRand(gapEnergy, E);
    if (qpE >= lowQPLimit*gapEnergy)
      qpEnergies.push_back(qpE);
    else
      EDep += qpE;

    if (E-qpE >= lowQPLimit*gapEnergy)
      qpEnergies.push_back(E-qpE);
    else
      EDep += E-qpE;
  }

  phonEnergies.clear();
  return EDep;
}

// Model the quasiparticles emitting phonons in the superconductor

G4double G4CMPProcessUtils::CalcPhononEnergies(G4double gapEnergy,
                                            G4double lowQPLimit,
                                            std::vector<G4double>& phonEnergies,
                                            std::vector<G4double>& qpEnergies) {
  // NOTE: Phonons with low energy will not be seen by the detector, so we
  // don't record those energies and just "lose" those phonons.
  // Have a reference in for loop b/c qp doesn't give all of its energy away.
  G4double EDep = 0.;
  std::vector<G4double> newQPEnergies;
  for (G4double& E: qpEnergies) {
    // NOTE: E mutates in PhononEnergyRand.
    G4double phonE = PhononEnergyRand(gapEnergy, E);
    if (phonE >= 2.0*gapEnergy)
      phonEnergies.push_back(phonE);
    if (E >= lowQPLimit*gapEnergy)
      newQPEnergies.push_back(E);
    else
      EDep += E;
  }

  qpEnergies.swap(newQPEnergies);
  return EDep;
}

// Calculate energies of phonon tracks that have reentered the crystal

void G4CMPProcessUtils::CalcReflectedPhononEnergies(
                                     G4MaterialPropertiesTable* prop,
                                     std::vector<G4double>& phonEnergies,
                                     std::vector<G4double>& reflectedEnergies) {
  // There is a 50% chance that a phonon is headed away from (toward) substrate
  std::vector<G4double> newPhonEnergies;
  for (G4double E : phonEnergies) {
    G4double frac = (G4UniformRand() < 0.5 ? 0.5 : 1.5);
    if (G4UniformRand() < CalcEscapeProbability(E, frac, prop))
      reflectedEnergies.push_back(E);
    else
      newPhonEnergies.push_back(E);
  }
  phonEnergies.swap(newPhonEnergies);
}

// Compute quasiparticle energy distribution from broken Cooper pair

G4double G4CMPProcessUtils::QPEnergyRand(G4double gapEnergy, G4double Energy) {
  // PDF is not integrable, so we can't do an inverse transform sampling.
  // Instead, we'll do a rejection method.
  //
  // PDF(E') = (E'*(Energy - E') + gapEnergy*gapEnergy)
  //           /
  //           sqrt((E'*E' - gapEnergy*gapEnergy) *
  //                ((Energy - E')*(Energy - E') - gapEnergy*gapEnergy));
  // The shape of the PDF is like a U, so the max values are at the endpoints:
  // E' = gapEnergy and E' = Energy - gapEnergy

  // Add buffer so first/last bins don't give zero denominator in pdfSum
  const G4double BUFF = 1000.;
  G4double xmin = gapEnergy + (Energy-2.*gapEnergy)/BUFF;
  G4double xmax = gapEnergy + (Energy-2.*gapEnergy)*(BUFF-1.)/BUFF;

  G4double ymax = (xmin*(Energy-xmin) + gapEnergy*gapEnergy)
                  /
                  sqrt((xmin*xmin - gapEnergy*gapEnergy) *
                       ((Energy-xmin)*(Energy-xmin) - gapEnergy*gapEnergy));

  G4double ytest = G4UniformRand()*ymax;
  G4double xtest = G4UniformRand()*(xmax-xmin) + xmin;
  while (ytest > (xtest*(Energy-xtest) + gapEnergy*gapEnergy)
                  /
                  sqrt((xtest*xtest - gapEnergy*gapEnergy) *
                       ((Energy-xtest)*(Energy-xtest) - gapEnergy*gapEnergy))) {
    ytest = G4UniformRand()*ymax;
    xtest = G4UniformRand()*(xmax-xmin) + xmin;
  }

  return xtest;
/*
  // PDF is not integrable, so we can't do an inverse transform sampling.
  // PDF is shaped like a capital U, so a rejection method would be slow.
  // Let's numerically calculate a CDF and throw a random to find E.

  const G4int BINS = 1000;
  G4double energyArr[BINS];
  G4double cdf[BINS];
  G4double pdfSum = 0.;
  for (size_t i = 0; i < BINS; ++i) {
    // Add 1 to i so first bin doesn't give zero denominator in pdfSum
    // Add 1 to BINS so last bin doesn't give zero denominator in pdfSum
    energyArr[i] = gapEnergy + (Energy-2.*gapEnergy) * (i+1.)/(BINS+1.);
    pdfSum += (energyArr[i]*(Energy - energyArr[i]) + gapEnergy*gapEnergy)
              /
              sqrt((energyArr[i]*energyArr[i] - gapEnergy*gapEnergy) *
                   ((Energy - energyArr[i])*(Energy - energyArr[i]) -
                    gapEnergy*gapEnergy));
    cdf[i] = pdfSum;
  }

  G4double u = G4UniformRand();

  size_t index = 0;
  for (; index < BINS; ++index) { //Combine normalization and search loops
    if (cdf[index]/pdfSum >= u)
      break;
  }

  return energyArr[index];
*/
}

// Compute phonon energy distribution from quasiparticle in superconductor

G4double
G4CMPProcessUtils::PhononEnergyRand(G4double gapEnergy, G4double& Energy) {
  // PDF is not integrable, so we can't do an inverse transform sampling.
  // Instead, we'll do a rejection method.
  //
  // PDF(E') = (E'*(Energy-E')*(Energy-E') * (E'-gapEnergy*gapEnergy/Energy))
  //           /
  //           sqrt((E'*E' - gapEnergy*gapEnergy);

  // Add buffer so first bin doesn't give zero denominator in pdfSum
  const G4double BUFF = 1000.;
  G4double xmin = gapEnergy + gapEnergy/BUFF;
  G4double xmax = Energy;

  G4double ymax = (xmin*(Energy-xmin)*(Energy-xmin) *
                    (xmin-gapEnergy*gapEnergy/Energy)) /
                  sqrt(xmin*xmin - gapEnergy*gapEnergy);

  G4double ytest = G4UniformRand()*ymax;
  G4double xtest = G4UniformRand()*(xmax-xmin) + xmin;
  while (ytest > (xtest*(Energy-xtest)*(Energy-xtest) *
                    (xtest-gapEnergy*gapEnergy/Energy)) /
                  sqrt(xtest*xtest - gapEnergy*gapEnergy)) {
    ytest = G4UniformRand()*ymax;
    xtest = G4UniformRand()*(xmax-xmin) + xmin;
  }

  G4double phononE = Energy - xtest;
  Energy = xtest;
  return phononE;
}


// Generate direction angle for phonons in Luke scattering

G4double G4CMPProcessUtils::MakePhononTheta(G4double k, G4double ks) const {
  G4double u = G4UniformRand();
  G4double v = ks/k;
  G4double base = (u-1) * (3*v - 3*v*v + v*v*v - 1);
  if (base < 0.0) return 0;
  
  G4double operand = v + pow(base, 1.0/3.0);   
  if (operand > 1.0) operand=1.0;
  
  return acos(operand);
}

// Compute energy of phonon in Luke Scattering

G4double G4CMPProcessUtils::MakePhononEnergy(G4double k, G4double ks,
					     G4double th_phonon) const {
  return 2.*(k*cos(th_phonon)-ks) * theLattice->GetSoundSpeed() * hbar_Planck;
}

// Compute direction angle for recoiling charge carrier

G4double G4CMPProcessUtils::MakeRecoilTheta(G4double k, G4double ks,
					    G4double th_phonon) const {
  if (th_phonon == 0.) return 0.;		// Avoid unnecessary work

  G4double kctks = k*cos(th_phonon) - ks;

  return acos( (k*k - 2*ks*kctks - 2*kctks*kctks)
	       / (k * sqrt(k*k - 4*ks*kctks)) );
}


// Construct new phonon or charge carrier track

G4Track* G4CMPProcessUtils::CreateTrack(G4ParticleDefinition* pd,
					const G4ThreeVector& waveVec,
					G4double energy) const {
  return CreateTrack(pd, waveVec, energy, currentTrack->GetPosition());
}

G4Track* G4CMPProcessUtils::CreateTrack(G4ParticleDefinition* pd,
					const G4ThreeVector& waveVec,
					G4double energy,
					const G4ThreeVector& pos) const {
  if (G4CMP::IsPhonon(pd)) {
    return CreatePhonon(G4PhononPolarization::Get(pd), waveVec, energy, pos);
  }

  if (G4CMP::IsChargeCarrier(pd)) {
    return CreateChargeCarrier(int(pd->GetPDGCharge()/eplus), ChooseValley(),
			       energy, waveVec, pos);
  }

  G4cerr << "WARNING: " << pd->GetParticleName() << " is not G4CMP" << G4endl;
  return new G4Track(new G4DynamicParticle(pd, waveVec, energy),
		     currentTrack->GetGlobalTime(), pos);
}


// Construct new phonon track with correct momentum, position, etc.

G4Track* G4CMPProcessUtils::CreatePhonon(G4int polarization,
					 const G4ThreeVector& waveVec,
					 G4double energy) const {
  return CreatePhonon(polarization,waveVec,energy,currentTrack->GetPosition());
}

G4Track* G4CMPProcessUtils::CreatePhonon(G4int polarization,
					 const G4ThreeVector& waveVec,
					 G4double energy,
           const G4ThreeVector& pos) const {
  if (polarization == G4PhononPolarization::UNKNOWN) {		// Choose value
    polarization = ChoosePhononPolarization();
  }

  G4ThreeVector vgroup = theLattice->MapKtoVDir(polarization, waveVec);
  if (std::fabs(vgroup.mag()-1.) > 0.01) {
    G4cerr << "WARNING: vgroup not a unit vector: " << vgroup
	   << " length " << vgroup.mag() << G4endl;
  }

  G4ParticleDefinition* thePhonon = G4PhononPolarization::Get(polarization);

  // Secondaries are (usually) created at the current track coordinates
  RotateToGlobalDirection(vgroup);
  G4ThreeVector secPos = AdjustSecondaryPosition(pos);

  G4Track* sec = new G4Track(new G4DynamicParticle(thePhonon, vgroup, energy),
                             currentTrack->GetGlobalTime(), secPos);

  // Store wavevector in auxiliary info for track
  AttachTrackInfo(sec)->SetPhononK(GetGlobalDirection(waveVec));

  sec->SetVelocity(theLattice->MapKtoV(polarization, waveVec));    
  sec->UseGivenVelocity(true);

  return sec;
}

// Generate random valley for charge carrier

G4int G4CMPProcessUtils::ChooseValley() const {
  return (G4int)(G4UniformRand()*theLattice->NumberOfValleys());
}

G4int G4CMPProcessUtils::ChangeValley(G4int valley) const {
  // generate random valley offset (up to N-1)
  G4int nv = theLattice->NumberOfValleys();
  G4int dv = (G4int)(G4UniformRand()*(nv-1))+1;

  // Apply offset to change input to new value
  return (valley+dv) % nv;
}


// Access electron propagation direction/index

G4int G4CMPProcessUtils::GetValleyIndex(const G4Track& track) const {
  return GetTrackInfo(track)->GetValleyIndex();
}

const G4RotationMatrix& 
G4CMPProcessUtils::GetValley(const G4Track& track) const {
  G4int iv = GetValleyIndex(track);
  return (iv>=0 ? theLattice->GetValley(iv) : G4RotationMatrix::IDENTITY);
}


// Construct new electron or hole track with correct conditions

G4Track* G4CMPProcessUtils::CreateChargeCarrier(G4int charge, G4int valley,
						const G4ThreeVector& p) const {
  return CreateChargeCarrier(charge, valley, p, currentTrack->GetPosition());
}

G4Track* 
G4CMPProcessUtils::CreateChargeCarrier(G4int charge, G4int valley,
				       G4double Ekin, 
				       const G4ThreeVector& dir,
				       const G4ThreeVector& pos) const {
  G4double carrierMass = 0.;
  if (charge==1)       carrierMass = theLattice->GetHoleMass();
  else if (charge==-1) carrierMass = theLattice->GetElectronMass();

  G4double carrierMom = std::sqrt(2.*Ekin*carrierMass);

  return CreateChargeCarrier(charge, valley, carrierMom*dir, pos);
}

G4Track* 
G4CMPProcessUtils::CreateChargeCarrier(G4int charge, G4int valley,
				       const G4ThreeVector& p,
				       const G4ThreeVector& pos) const {
  if (charge != 1 && charge != -1) {
    G4cerr << "ERROR:  CreateChargeCarrier invalid charge " << charge << G4endl;
    return 0;
  }

  G4ParticleDefinition* theCarrier = 0;
  G4double carrierMass=0., carrierEnergy=0.;

  G4ThreeVector v_unit;
  if (charge==1) {
    theCarrier    = G4CMPDriftHole::Definition();
    carrierMass   = theLattice->GetHoleMass();
    carrierEnergy = 0.5 * p.mag2() / carrierMass;	// Non-relativistic
    v_unit = p.unit();
  } else {
    theCarrier    = G4CMPDriftElectron::Definition();
    carrierMass   = theLattice->GetElectronMass();
    G4ThreeVector p_local = GetLocalDirection(p);
    G4ThreeVector v_local = theLattice->MapPtoV_el(valley, p_local);
    RotateToGlobalDirection(v_local);
    carrierEnergy = 0.5 * carrierMass * v_local.mag2();// Non-relativistic
    v_unit = v_local.unit();
  }

  G4DynamicParticle* secDP =
    new G4DynamicParticle(theCarrier, v_unit, carrierEnergy, carrierMass);

  G4ThreeVector secPos = AdjustSecondaryPosition(pos);
  G4Track* sec = new G4Track(secDP, currentTrack->GetGlobalTime(), secPos);

  // Store wavevector in auxiliary info for track
  AttachTrackInfo(sec)->SetValleyIndex(valley);

  return sec;
}

G4ThreeVector G4CMPProcessUtils::AdjustSecondaryPosition(G4ThreeVector pos) const {
  // Take a copy because we would've had to make a copy at some point anyway.
  // If the step is near a boundary, create the secondary in the initial volume
  G4Navigator* nav = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  // Safety is the distance to a boundary
  G4double safety = nav->ComputeSafety(pos);
  // Tolerance is the error in deciding which volume a track is in
  G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  // If the distance to an edge is within error, we might accidentally get placed
  // in the next volume. Instead, let's scoot a bit away from the edge.
  if (safety <= kCarTolerance) {
    G4ThreeVector norm = currentVolume->GetLogicalVolume()->GetSolid()->SurfaceNormal(GetLocalPosition(pos));
    RotateToGlobalDirection(norm);
    pos += (safety - kCarTolerance * (1.001)) * norm;
  }

  return pos;
}

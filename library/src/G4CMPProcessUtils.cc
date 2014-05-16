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

#include "G4CMPProcessUtils.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4CMPValleyTrackMap.hh"
#include "G4AffineTransform.hh"
#include "G4DynamicParticle.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhononLong.hh"
#include "G4PhononPolarization.hh"
#include "G4PhononTrackMap.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononTransSlow.hh"
#include "G4PhysicalConstants.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4VTouchable.hh"
#include "Randomize.hh"


// Constructor and destructor

G4CMPProcessUtils::G4CMPProcessUtils()
  : theLattice(0), trackKmap(G4PhononTrackMap::GetInstance()),
    trackVmap(G4CMPValleyTrackMap::GetInstance()), currentTrack(0) {;}

G4CMPProcessUtils::~G4CMPProcessUtils() {;}


// Initialization for current track

void G4CMPProcessUtils::LoadDataForTrack(const G4Track* track) {
  currentTrack = track;

  // WARNING!  This assumes track starts and ends in one single volume!
  FindLattice(track->GetVolume());
  SetTransforms(track->GetTouchable());

  // Register track in either phonon or charge-carrier map
  const G4ParticleDefinition* pd = track->GetParticleDefinition();

  if (pd == G4PhononLong::Definition() ||
      pd == G4PhononTransFast::Definition() ||
      pd == G4PhononTransSlow::Definition()) {
    if (!trackKmap->Find(track)) 
      // FIXME:  THE WAVEVECTOR SHOULD BE COMPUTED BY INVERTING THE K/V MAP
      trackKmap->SetK(track, track->GetMomentumDirection());
  }

  if (pd == G4CMPDriftElectron::Definition() ||
      pd == G4CMPDriftHole::Definition()) {
    if (!trackVmap->Find(track))
      // FIXME:  HOW DO WE CONVERT THE MOMENTUM TO AN INITIAL VALLEY?
      trackVmap->SetValley(track, ChooseValley());
  }
}

// Fetch lattice for current track, use in subsequent steps

void G4CMPProcessUtils::FindLattice(const G4VPhysicalVolume* volume) {
  G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();
  theLattice = LM->GetLattice(volume);
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
  SetTransforms(0);
  trackKmap->RemoveTrack(currentTrack);
  trackVmap->RemoveTrack(currentTrack);
  currentTrack = 0;
  theLattice = 0;
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
  return GetLocalDirection(track.GetMomentum());
}

void G4CMPProcessUtils::GetLocalMomentum(const G4Track& track, 
					 G4double mom[3]) const {
  G4ThreeVector tmom = GetLocalMomentum(track);
  mom[0] = tmom.x();
  mom[1] = tmom.y();
  mom[2] = tmom.z();
}


// Return particle type for currently active track [set in LoadDataForTrack()]

const G4ParticleDefinition* G4CMPProcessUtils::GetCurrentParticle() const {
  return (currentTrack ? currentTrack->GetParticleDefinition() : 0);
}


// Access phonon particle-type/polarization indices

G4int G4CMPProcessUtils::GetPolarization(const G4Track& track) const {
  return G4PhononPolarization::Get(track.GetParticleDefinition());
}


// Generate random polarization from density of states

G4int G4CMPProcessUtils::ChoosePolarization(G4double Ldos, G4double STdos,
					    G4double FTdos) const {
  G4double norm = Ldos + STdos + FTdos;
  G4double cProbST = STdos/norm;
  G4double cProbFT = FTdos/norm + cProbST;

  // NOTE:  Order of selection done to match previous random sequences
  G4double modeMixer = G4UniformRand();
  if (modeMixer<cProbST) return G4PhononPolarization::TransSlow;
  if (modeMixer<cProbFT) return G4PhononPolarization::TransFast;
  return G4PhononPolarization::Long;
}

G4int G4CMPProcessUtils::ChoosePolarization() const {
  return ChoosePolarization(theLattice->GetLDOS(), theLattice->GetSTDOS(),
			    theLattice->GetFTDOS());
}


// Construct new phonon track with correct momentum, position, etc.

G4Track* G4CMPProcessUtils::CreatePhonon(G4int polarization,
					 const G4ThreeVector& waveVec,
					 G4double energy) const {
  if (polarization == G4PhononPolarization::UNKNOWN) {		// Choose value
    polarization = ChoosePolarization();
  }

  G4ThreeVector vgroup = theLattice->MapKtoVDir(polarization, waveVec);
  if (std::fabs(vgroup.mag()-1.) > 0.01) {
    G4cerr << "WARNING: vgroup not a unit vector: " << vgroup
	   << " length " << vgroup.mag() << G4endl;
  }

  G4ParticleDefinition* thePhonon = G4PhononPolarization::Get(polarization);

  // Secondaries are created at the current track coordinates
  RotateToGlobalDirection(vgroup);
  G4Track* sec = new G4Track(new G4DynamicParticle(thePhonon, vgroup, energy),
			     currentTrack->GetGlobalTime(),
			     currentTrack->GetPosition());

  // Store wavevector in lookup table for future tracking
  trackKmap->SetK(sec, GetGlobalDirection(waveVec));

  sec->SetVelocity(theLattice->MapKtoV(polarization, waveVec));    
  sec->UseGivenVelocity(true);

  return sec;
}


// Generate random valley for charge carrier

G4int G4CMPProcessUtils::ChooseValley() const {
  return (G4int)(G4UniformRand()*theLattice->NumberOfValleys());  
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
  if (th_phonon == 0.) return 0.;		// Avoid unnecessary work

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


// Access electron propagation direction/index

G4int G4CMPProcessUtils::GetValleyIndex(const G4Track& track) const {
  return trackVmap->GetValley(&track);
}

const G4RotationMatrix& 
G4CMPProcessUtils::GetValley(const G4Track& track) const {
  G4int iv = GetValleyIndex(track);
  return (iv>=0 ? theLattice->GetValley(iv) : G4RotationMatrix::IDENTITY);
}


// Construct new electron or hole track with correct conditions

G4Track* G4CMPProcessUtils::CreateChargeCarrier(G4int charge, G4int valley,
						const G4ThreeVector& p) const {
  if (charge != 1 && charge != -1) {
    G4cerr << "ERROR:  CreateChargeCarrier invalid charge " << charge << G4endl;
    return 0;
  }

  G4ParticleDefinition* theCarrier = 0;
  G4double carrierMass=0., carrierEnergy=0.;
  if (charge==1) {
    theCarrier    = G4CMPDriftHole::Definition();
    carrierMass   = theLattice->GetHoleMass();
    carrierEnergy = 0.5 * p.mag2() / carrierMass;	// Non-relativistic
  } else {
    theCarrier    = G4CMPDriftElectron::Definition();
    G4ThreeVector p_local = GetLocalDirection(p);
    carrierMass   = theLattice->GetElectronEffectiveMass(valley, p_local);
    carrierEnergy = theLattice->MapPtoEkin(valley, p_local);
  }

  G4DynamicParticle* secDP = new G4DynamicParticle(theCarrier,p,carrierEnergy);
  secDP->SetMass(carrierMass);

  G4Track* sec = new G4Track(secDP, currentTrack->GetGlobalTime(),
			     currentTrack->GetPosition());

  // Store valley index in lookup table for future tracking
  trackVmap->SetValley(sec, valley);

  return sec;
}

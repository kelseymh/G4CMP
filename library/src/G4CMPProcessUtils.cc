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

  // Fetch lattice for current track once, use in subsequent steps
  // WARNING!  This assumes track starts and ends in one single volume!
  G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();
  theLattice = LM->GetLattice(track->GetVolume());

  // Get coordinate transformations for current volume
  // WARNING!  This assumes track starts and ends in one single volume!
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

void G4CMPProcessUtils::SetTransforms(const G4VTouchable* touchable) {
  if (!touchable) {			// Null pointer defaults to identity
    fLocalToGlobal = fGlobalToLocal = G4AffineTransform();
    return;
  }

  fLocalToGlobal = G4AffineTransform(touchable->GetRotation(),
				     touchable->GetTranslation());
  fGlobalToLocal = fLocalToGlobal.Inverse();
}

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


// Construct new phonon track with correct momentum, position, etc.

G4Track* G4CMPProcessUtils::CreatePhonon(G4int polarization,
					 const G4ThreeVector& waveVec,
					 G4double energy) const {
  G4ThreeVector vgroup = theLattice->MapKtoVDir(polarization, waveVec);
  if (std::fabs(vgroup.mag()-1.) > 0.01) {
    G4cout << "WARNING: vgroup not a unit vector: " << vgroup << G4endl;
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
						const G4ThreeVector& p,
						G4double energy) const {
  if (charge != 1 && charge != -1) {
    G4cerr << "ERROR:  CreateChargeCarrier invalid charge " << charge << G4endl;
    return 0;
  }

  G4ParticleDefinition* theCarrier = 0;
  if (charge==1) theCarrier = G4CMPDriftHole::Definition();
  else theCarrier = G4CMPDriftElectron::Definition();

  G4Track* sec = new G4Track(new G4DynamicParticle(theCarrier, p, energy),
			     currentTrack->GetGlobalTime(),
			     currentTrack->GetPosition());

  // Store valley index in lookup table for future tracking
  trackVmap->SetValley(sec, valley);

  // FIXME:  Should momentum be aligned with valley here?

  return sec;
}

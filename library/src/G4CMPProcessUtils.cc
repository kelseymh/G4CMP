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

#include "G4CMPProcessUtils.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4Track.hh"
#include "G4PhononTrackMap.hh"
#include "G4CMPValleyTrackMap.hh"
#include "G4PhononPolarization.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4CMPDriftHole.hh"
#include "G4CMPDriftElectron.hh"


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
}


// Access phonon particle-type/polarization indices

G4int G4CMPProcessUtils::GetPolarization(const G4Track& track) const {
  return G4PhononPolarization::Get(track.GetParticleDefinition());
}


// Construct new phonon track with correct momentum, position, etc.

G4Track* G4CMPProcessUtils::CreatePhonon(G4int polarization,
					 const G4ThreeVector& waveVec,
					 G4double energy) const {
  G4ThreeVector vgroup = theLattice->MapKtoVDir(polarization, waveVec);
  vgroup = theLattice->RotateToGlobal(vgroup);
  if (std::fabs(vgroup.mag()-1.) > 0.01) {
    G4cout << "WARNING: vgroup not a unit vector: " << vgroup << G4endl;
  }

  G4ParticleDefinition* thePhonon = G4PhononPolarization::Get(polarization);

  // Secondaries are created at the current track coordinates
  G4Track* sec = new G4Track(new G4DynamicParticle(thePhonon, vgroup, energy),
			     currentTrack->GetGlobalTime(),
			     currentTrack->GetPosition());

  // Store wavevector in lookup table for future tracking
  trackKmap->SetK(sec, theLattice->RotateToGlobal(waveVec));

  sec->SetVelocity(theLattice->MapKtoV(polarization, waveVec));    
  sec->UseGivenVelocity(true);

  return sec;
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
  if (charge != 0 && charge != -1) {
    G4cerr << "ERROR:  CreateChargeCarrier invalid charge " << charge << G4endl;
    return 0;
  }

  G4ParticleDefinition* theCarrier = 0;
  if (charge==0) theCarrier = G4CMPDriftHole::Definition();
  else theCarrier = G4CMPDriftElectron::Definition();

  G4Track* sec = new G4Track(new G4DynamicParticle(theCarrier, p, energy),
			     currentTrack->GetGlobalTime(),
			     currentTrack->GetPosition());

  // Store valley index in lookup table for future tracking
  trackVmap->SetValley(sec, valley);

  // FIXME:  Should momentum be aligned with valley here?

  return sec;
}

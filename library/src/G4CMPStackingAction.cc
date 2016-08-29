/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPStackingAction.cc
/// \brief Implementation of the G4CMPStackingAction class
///     This stacking action is necessary to ensure that velocity and 
///     propagation direction are set properly for phonons created with
///     G4ParticleGun, and to ensure that the initial lattice valley
///	is set properly for created drifting electrons.
//
// $Id$
//
// 20140411 Set charge carrier masses appropriately for material
// 20141216 Set velocity for electrons
// 20150109 Protect velocity flag with compiler flag
// 20160625 Process _all_ tracks to ensure they're in correct volumes
// 20160829  Drop G4CMP_SET_ELECTRON_MASS code blocks; not physical

#include "G4CMPStackingAction.hh"
#include "G4CMPTrackInformation.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4PhononLong.hh"
#include "G4PhononPolarization.hh"
#include "G4PhononTrackMap.hh"
#include "G4CMPDriftHole.hh"
#include "G4CMPDriftElectron.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononTransSlow.hh"
#include "G4PhysicalConstants.hh"
#include "G4RandomDirection.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "Randomize.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4CMPStackingAction::G4CMPStackingAction()
  : G4UserStackingAction(), G4CMPProcessUtils() {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4CMPStackingAction::~G4CMPStackingAction() {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ClassificationOfNewTrack 
G4CMPStackingAction::ClassifyNewTrack(const G4Track* aTrack) {
  G4ClassificationOfNewTrack classification = fUrgent;

  // Configure utility functions for current track
  FindLattice(aTrack->GetVolume());
  SetTransforms(aTrack->GetTouchable());

  // If phonon or charge carrier is not in a lattice-enabled volume, kill it
  if ((IsPhonon(aTrack) || IsChargeCarrier(aTrack)) && !theLattice)
    return fKill;
  
  // Non-initial tracks should not be touched
  if (aTrack->GetParentID() != 0) return classification;

  // Attach auxiliary info to new track
  AttachTrackInfo(aTrack);

  // Fill kinematic data for new track (secondaries will have this done)
  if (IsPhonon(aTrack)) {
    SetPhononWaveVector(aTrack);
    SetPhononVelocity(aTrack);
  }

  if (IsChargeCarrier(aTrack)) {
    SetChargeCarrierValley(aTrack);
    SetChargeCarrierMass(aTrack);
  }

  ReleaseTrack();

  return classification; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Generate random wave vector for initial phonon (ought to invert momentum)

void G4CMPStackingAction::SetPhononWaveVector(const G4Track* aTrack) const {
  //Compute random wave-vector (override whatever ParticleGun did)
  G4ThreeVector ranK = G4RandomDirection();
  
  //Store wave-vector as track information
  AttachTrackInfo(aTrack)->SetPhononK(ranK);
}

// Set velocity of phonon track appropriately for material

void G4CMPStackingAction::SetPhononVelocity(const G4Track* aTrack) const {
  // Get wavevector associated with track
  G4ThreeVector K = GetTrackInfo(aTrack)->GetPhononK();
  G4int pol = GetPolarization(aTrack);

  //Compute direction of propagation from wave vector
  G4ThreeVector momentumDir = theLattice->MapKtoVDir(pol, K);

  if (momentumDir.mag() < 0.9) {
    G4cerr << " track mode " << pol << " K " << K << G4endl;
    G4Exception("G4CMPStackingAction::SetPhononVelocity", "Lattice010",
		FatalException, "KtoVDir failed to return unit vector");
    return;
  }

  //Compute true velocity of propagation
  G4double velocity = theLattice->MapKtoV(pol, K);
  
  // Cast to non-const pointer so we can adjust non-standard kinematics
  G4Track* theTrack = const_cast<G4Track*>(aTrack);

  theTrack->SetMomentumDirection(momentumDir);
  theTrack->SetVelocity(velocity);
  theTrack->UseGivenVelocity(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Set current Brillouin valley (random) for electrons in material

void G4CMPStackingAction::SetChargeCarrierValley(const G4Track* aTrack) const {
  if (aTrack->GetDefinition() != G4CMPDriftElectron::Definition()) return;

  if (GetValleyIndex(aTrack) < 0) {
    int valley = (G4int)(G4UniformRand()*theLattice->NumberOfValleys());
    AttachTrackInfo(aTrack)->SetValleyIndex(valley);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Set dynamical mass of charge carrier to scalar value for material

void G4CMPStackingAction::SetChargeCarrierMass(const G4Track* aTrack) const {
  // Get effective mass for charge carrier
  G4double mass = aTrack->GetDefinition()->GetPDGMass();

  if (IsHole(aTrack))     mass = theLattice->GetHoleMass();
  if (IsElectron(aTrack)) mass = theLattice->GetElectronMass();	// H-V scalar

  // Cast to non-const pointer so we can change the effective mass
  G4DynamicParticle* dynp =
    const_cast<G4DynamicParticle*>(aTrack->GetDynamicParticle());

  dynp->SetMass(mass*c_squared);	// Converts to Geant4 [M]=[E] units
}

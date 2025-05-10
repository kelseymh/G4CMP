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
// 20160829 Drop G4CMP_SET_ELECTRON_MASS code blocks; not physical
// 20170620 Drop obsolete SetTransforms() call
// 20170624 Clean up track initialization
// 20170928 Replace "polarization" with "mode"
// 20211001 Remove electron energy adjustment; set mass instead.
//		Assign electron valley nearest to momentum direction.
// 20230702 I. Ataee -- Corrections to effective mass calculations
//		for new charge carriers.
// 20240122 G4CMP-446 -- SetPhononVelocity() should use global-to-local
//		transform for k vector and Vg.
// 20250508 N. Tenpas -- Add coordinate transforms in SetPhononVelocity.

#include "G4CMPStackingAction.hh"

#include "G4CMPDriftHole.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftTrackInfo.hh"
#include "G4CMPPhononTrackInfo.hh"
#include "G4CMPTrackUtils.hh"
#include "G4CMPUtils.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4PhononLong.hh"
#include "G4PhononPolarization.hh"
#include "G4PhononTrackMap.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononTransSlow.hh"
#include "G4PhysicalConstants.hh"
#include "G4RandomDirection.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
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

  // Don't do anything to ordinary G4 tracks
  if (!G4CMP::IsPhonon(aTrack) && !G4CMP::IsChargeCarrier(aTrack))
    return classification;

  // Configure utility functions for current track (do NOT use LoadDataForTrack)
  SetCurrentTrack(aTrack);
  SetLattice(aTrack);

  // If phonon or charge carrier is not in a lattice-enabled volume, kill track immediately
  if ((IsPhonon() || IsChargeCarrier()) && !theLattice) {
    ReleaseTrack();
    return fKill;
  }

  // Attach appropriate container to store additional kinematics if needed
  if (!G4CMP::HasTrackInfo(aTrack)) {
    G4CMP::AttachTrackInfo(aTrack);

    // Fill kinematic data for new track (secondaries will have this done)
    if (IsPhonon()) SetPhononVelocity(aTrack);

    if (IsChargeCarrier()) {
      AssignNearestValley(aTrack);
      SetChargeCarrierMass(aTrack);
    }
  }

  ReleaseTrack();

  return classification; 
}

// Set velocity of phonon track appropriately for material

void G4CMPStackingAction::SetPhononVelocity(const G4Track* aTrack) const {
  // Get wavevector associated with track
  G4ThreeVector k = G4CMP::GetTrackInfo<G4CMPPhononTrackInfo>(*aTrack)->k();
  G4int mode = GetPolarization(aTrack);

  // Compute direction of propagation from wave vector
  // Geant4 thinks that momentum and velocity point in same direction,
  // momentumDir here actually means velocity direction.

  RotateToLocalDirection(k);	// G4LatticePhysical expects local-frame vector
  G4ThreeVector momentumDir = theLattice->MapKtoVDir(mode, k);
  RotateToGlobalDirection(momentumDir);	      // and returns local-frame vector

  if (momentumDir.mag() < 0.9) {
    G4cerr << " track mode " << mode << " k " << k << G4endl;
    G4Exception("G4CMPStackingAction::SetPhononVelocity", "Lattice010",
		FatalException, "KtoVDir failed to return unit vector");
    return;
  }

  // Compute true velocity of propagation
  G4double velocity = theLattice->MapKtoV(mode, k);
  
  // Cast to non-const pointer so we can adjust non-standard kinematics
  G4Track* theTrack = const_cast<G4Track*>(aTrack);

  theTrack->SetMomentumDirection(momentumDir);
  theTrack->SetVelocity(velocity);
  theTrack->UseGivenVelocity(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Assign electron to valley nearest to momentum direction

void G4CMPStackingAction::AssignNearestValley(const G4Track* aTrack) const {
  G4int valley = FindNearestValley(aTrack);

  G4CMP::GetTrackInfo<G4CMPDriftTrackInfo>(*aTrack)->SetValleyIndex(valley);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Set dynamical mass of charge carrier to scalar value for material

void G4CMPStackingAction::SetChargeCarrierMass(const G4Track* aTrack) const {
  if (!G4CMP::IsChargeCarrier(aTrack)) return;

  G4int iv = GetCurrentValley();
  G4ThreeVector pdir = aTrack->GetMomentumDirection();
  G4double ekin = aTrack->GetKineticEnergy();
  G4ThreeVector p = theLattice->MapEkintoP(iv,GetLocalDirection(pdir),ekin);

  G4double mass = 
    G4CMP::IsHole(aTrack) ? theLattice->GetHoleMass() :
    G4CMP::IsElectron(aTrack) ? 
    theLattice->GetElectronEffectiveMass(iv,p) :
    aTrack->GetDynamicParticle()->GetMass()/c_squared;

  // Cast to non-const pointer so we can change the effective mass
  G4DynamicParticle* dynp =
    const_cast<G4DynamicParticle*>(aTrack->GetDynamicParticle());

  dynp->SetMass(mass*c_squared);	// Converts to Geant4 [M]=[E] units
}

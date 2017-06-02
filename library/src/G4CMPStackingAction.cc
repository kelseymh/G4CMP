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
  if ((G4CMP::IsPhonon(aTrack) || G4CMP::IsChargeCarrier(aTrack))
      && !theLattice)
    return fKill;
  
  // Non-initial tracks should not be touched
  if (aTrack->GetParentID() != 0) return classification;

  // Fill kinematic data for new track (secondaries will have this done)
  if (G4CMP::IsPhonon(aTrack)) {
    auto trackInfo = new G4CMPPhononTrackInfo(theLattice, G4RandomDirection());
    G4CMP::AttachTrackInfo(*aTrack, trackInfo);
    SetPhononVelocity(aTrack);
  }

  if (G4CMP::IsChargeCarrier(aTrack)) {
    SetChargeCarrierMass(aTrack);
    if (G4CMP::IsElectron(aTrack)) {
      auto trackInfo = new G4CMPDriftTrackInfo(theLattice, G4CMP::ChooseValley(theLattice));
      G4CMP::AttachTrackInfo(*aTrack, trackInfo);
      SetElectronEnergy(aTrack);
    } else { // IsHole
      auto trackInfo = new G4CMPDriftTrackInfo(theLattice, -1);
      G4CMP::AttachTrackInfo(*aTrack, trackInfo);
    }
  }

  ReleaseTrack();

  return classification; 
}

// Set velocity of phonon track appropriately for material

void G4CMPStackingAction::SetPhononVelocity(const G4Track* aTrack) const {
  // Get wavevector associated with track
  G4ThreeVector k = G4CMP::GetTrackInfo<G4CMPPhononTrackInfo>(*aTrack)->k();
  G4int pol = GetPolarization(aTrack);

  // Compute direction of propagation from wave vector
  // Geant4 thinks that momentum and velocity point in same direction,
  // momentumDir here actually means velocity direction.
  G4ThreeVector momentumDir = theLattice->MapKtoVDir(pol, k);

  if (momentumDir.mag() < 0.9) {
    G4cerr << " track mode " << pol << " k " << k << G4endl;
    G4Exception("G4CMPStackingAction::SetPhononVelocity", "Lattice010",
		FatalException, "KtoVDir failed to return unit vector");
    return;
  }

  //Compute true velocity of propagation
  G4double velocity = theLattice->MapKtoV(pol, k);
  
  // Cast to non-const pointer so we can adjust non-standard kinematics
  G4Track* theTrack = const_cast<G4Track*>(aTrack);

  theTrack->SetMomentumDirection(momentumDir);
  theTrack->SetVelocity(velocity);
  theTrack->UseGivenVelocity(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Set dynamical mass of charge carrier to scalar value for material

void G4CMPStackingAction::SetChargeCarrierMass(const G4Track* aTrack) const {
  // Get effective mass for charge carrier
  G4double mass = aTrack->GetDefinition()->GetPDGMass();

  if (G4CMP::IsHole(aTrack))     mass = theLattice->GetHoleMass();
  if (G4CMP::IsElectron(aTrack)) mass = theLattice->GetElectronMass();	// H-V scalar

  // Cast to non-const pointer so we can change the effective mass
  G4DynamicParticle* dynp =
    const_cast<G4DynamicParticle*>(aTrack->GetDynamicParticle());

  dynp->SetMass(mass*c_squared);	// Converts to Geant4 [M]=[E] units
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Set G4Track energy to correctly calculate velocity

void G4CMPStackingAction::SetElectronEnergy(const G4Track* aTrack) const {
  G4int valley = G4CMP::GetTrackInfo<G4CMPDriftTrackInfo>(*aTrack)->ValleyIndex();
  G4double E = aTrack->GetKineticEnergy();
  G4double kmag_HV = std::sqrt(2. * E * theLattice->GetElectronMass()) /
                     hbar_Planck;
  G4ThreeVector kdir_HV = theLattice->MapV_elToK_HV(valley,
                                                    aTrack->GetMomentumDirection());

  G4ThreeVector p = theLattice->MapK_HVtoP(valley, kmag_HV * kdir_HV.unit());
  G4ThreeVector vTrue = theLattice->MapPtoV_el(valley, p);
  // Set fake E that yeilds correct v
  (const_cast<G4Track*>(aTrack))->SetKineticEnergy(0.5 *
                                                   aTrack->GetDynamicParticle()->GetMass() /
                                                   c_squared *
                                                   vTrue.mag2());
}

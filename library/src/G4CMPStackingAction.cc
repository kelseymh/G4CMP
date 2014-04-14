//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
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

#include "G4CMPStackingAction.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4PhononLong.hh"
#include "G4PhononPolarization.hh"
#include "G4PhononTrackMap.hh"
#include "G4CMPDriftHole.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPValleyTrackMap.hh"
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

G4CMPStackingAction::G4CMPStackingAction() {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4CMPStackingAction::~G4CMPStackingAction() {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ClassificationOfNewTrack 
G4CMPStackingAction::ClassifyNewTrack(const G4Track* aTrack) {
  G4ClassificationOfNewTrack classification = fUrgent;

  // Non-initial tracks should not be touched
  if (aTrack->GetParentID() != 0) return classification;

  G4ParticleDefinition* pd = aTrack->GetDefinition();

  // Obtain LatticeManager and lattice for current volume
  // WARNING!  This assumes track starts and ends in one single volume!
  G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();
  G4LatticePhysical* theLattice = LM->GetLattice(aTrack->GetVolume());

  if (pd == G4PhononLong::Definition() ||
      pd == G4PhononTransFast::Definition() ||
      pd == G4PhononTransSlow::Definition()) {
    G4int pol = G4PhononPolarization::Get(aTrack->GetDefinition());
    
    //Compute random wave-vector (override whatever ParticleGun did)
    G4ThreeVector Ran = G4RandomDirection();
    
    //Store wave-vector as track information
    G4PhononTrackMap* theKmap = G4PhononTrackMap::GetPhononTrackMap();
    theKmap->SetK(aTrack, Ran);
    
    //Compute direction of propagation from wave vector
    G4ThreeVector momentumDir = theLattice->MapKtoVDir(pol, Ran);
    
    //Compute true velocity of propagation
    G4double velocity = theLattice->MapKtoV(pol, Ran);
    
    //cast to non-const pointer so we can set the velocity
    G4Track* theTrack = const_cast<G4Track*>(aTrack);

    theTrack->SetMomentumDirection(momentumDir);
    theTrack->SetVelocity(velocity);
    theTrack->UseGivenVelocity(true);
  }

  if (pd == G4CMPDriftHole::Definition() ||
      pd == G4CMPDriftElectron::Definition()) {
    SetChargeCarrierMass(aTrack);
    SetChargeCarrierValley(aTrack);
  }

  return classification; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Set dynamical mass of charge carrier to scalar value for material

void G4CMPStackingAction::SetChargeCarrierMass(const G4Track* aTrack) const {
  G4ParticleDefinition* pd = aTrack->GetDefinition();
  if (pd != G4CMPDriftHole::Definition() &&
      pd != G4CMPDriftElectron::Definition()) return;

  // Obtain LatticeManager and lattice for current volume
  // WARNING!  This assumes track starts and ends in one single volume!
  G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();
  G4LatticePhysical* theLattice = LM->GetLattice(aTrack->GetVolume());

  // Cast to non-const pointer so we can change the effective mass
  G4DynamicParticle* dynp =
    const_cast<G4DynamicParticle*>(aTrack->GetDynamicParticle());

  G4double mass = (pd == G4CMPDriftHole::Definition()
		   ? theLattice->GetHoleMass() : theLattice->GetElectronMass());
  dynp->SetMass(mass*c_squared);	// Converts to Geant4 [M]=[E] units
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Set current Brillouin valley (random) for electrons in material

void G4CMPStackingAction::SetChargeCarrierValley(const G4Track* aTrack) const {
  if (aTrack->GetDefinition() != G4CMPDriftElectron::Definition()) return;

  // Obtain LatticeManager and lattice for current volume
  // WARNING!  This assumes track starts and ends in one single volume!
  G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();
  G4LatticePhysical* theLattice = LM->GetLattice(aTrack->GetVolume());

  int valley = (G4int)(G4UniformRand()*theLattice->NumberOfValleys());
  G4CMPValleyTrackMap* theIVmap = G4CMPValleyTrackMap::GetValleyTrackMap();
  theIVmap->SetValley(aTrack, valley);
}

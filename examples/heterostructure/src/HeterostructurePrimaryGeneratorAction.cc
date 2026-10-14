/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file HeterostructurePrimaryGeneratorAction.cc 
/// \brief Implementation of the HeterostructurePrimaryGeneratorAction class
/// Constructs heterostructure toy geometry with germanium on silicon.

// 20260815  Selby Q. Dang (Stanford/SLAC)

#include "HeterostructurePrimaryGeneratorAction.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4Event.hh"
#include "G4Geantino.hh"
#include "G4ParticleGun.hh"
#include "G4RandomDirection.hh"
#include "G4SystemOfUnits.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

HeterostructurePrimaryGeneratorAction::HeterostructurePrimaryGeneratorAction() { 
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);   

  // default particle kinematics 
  fParticleGun->SetParticleDefinition(G4Geantino::Definition());
  fParticleGun->SetParticleMomentumDirection(G4RandomDirection());
  fParticleGun->SetParticlePosition(G4ThreeVector(0.0,0.0,90.*nm));
  fParticleGun->SetParticleEnergy(1e-6*eV);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

HeterostructurePrimaryGeneratorAction::~HeterostructurePrimaryGeneratorAction() {
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HeterostructurePrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
  if (fParticleGun->GetParticleDefinition() == G4Geantino::Definition()) {
    // Generate hole
    fParticleGun->SetParticleDefinition(G4CMPDriftHole::Definition());
    fParticleGun->SetParticleMomentumDirection(G4RandomDirection());
    fParticleGun->GeneratePrimaryVertex(anEvent);
    
    // Generate electron with different random direction
    fParticleGun->SetParticleDefinition(G4CMPDriftElectron::Definition());
    fParticleGun->SetParticleMomentumDirection(G4RandomDirection());
    fParticleGun->GeneratePrimaryVertex(anEvent);

    // Restore "not set" condition for next event
    fParticleGun->SetParticleDefinition(G4Geantino::Definition());
  } else {
    fParticleGun->GeneratePrimaryVertex(anEvent);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



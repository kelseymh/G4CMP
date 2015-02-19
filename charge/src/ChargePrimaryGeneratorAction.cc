// $Id$
//
// Generator uses G4ParticleGun, producing one electron and one hole per
// event by default.  User may change the number of particles per event
// via macro command |/gun/number|, the starting position, or the energy.

#include "ChargePrimaryGeneratorAction.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4Geantino.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"


ChargePrimaryGeneratorAction::ChargePrimaryGeneratorAction() {
  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);
  
  // default particle kinematics -- user may specify individual particle
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1));
  particleGun->SetParticlePosition(G4ThreeVector(0.0,0.0,0.0));
  particleGun->SetParticleEnergy(1e-13*eV);
}

ChargePrimaryGeneratorAction::~ChargePrimaryGeneratorAction() {
  delete particleGun;
}

void ChargePrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
  // If user did not set particle explicitly, do e/h pairs
  if (particleGun->GetParticleDefinition() == G4Geantino::Definition()) {
    particleGun->SetParticleDefinition(G4CMPDriftHole::Definition());
    particleGun->GeneratePrimaryVertex(anEvent);
    particleGun->SetParticleDefinition(G4CMPDriftElectron::Definition());
    particleGun->GeneratePrimaryVertex(anEvent);

    // Restore "not set" condition for next event
    particleGun->SetParticleDefinition(G4Geantino::Definition());
  } else {
    particleGun->GeneratePrimaryVertex(anEvent);
  }
}


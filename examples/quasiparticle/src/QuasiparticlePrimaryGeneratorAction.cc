/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file exoticphysics/quasiparticle/src/QuasiparticlePrimaryGeneratorAction.cc
/// \brief Implementation of the QuasiparticlePrimaryGeneratorAction class
//
// $Id: e75f788b103aef810361fad30f75077829192c13 $
//
// 20140519  Allow the user to specify phonon type by name in macro; if
//	     "geantino" is set, use random generator to select.

#include "QuasiparticlePrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4Geantino.hh"
#include "G4ParticleGun.hh"
#include "G4RandomDirection.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononTransSlow.hh"
#include "G4PhononLong.hh"
#include "G4SystemOfUnits.hh"
#include "G4GeneralParticleSource.hh"

using namespace std;

QuasiparticlePrimaryGeneratorAction::QuasiparticlePrimaryGeneratorAction() { 
  //G4int n_particle = 1;
  //fParticleGun  = new G4ParticleGun(n_particle);   

  // default particle kinematics ("geantino" triggers random phonon choice)
  //fParticleGun->SetParticleDefinition(G4Geantino::Definition());
  //fParticleGun->SetParticleMomentumDirection(G4RandomDirection());
  //fParticleGun->SetParticlePosition(G4ThreeVector(0.0,0.0,0.0));
  //fParticleGun->SetParticleEnergy(0.0075*eV);

  fParticleGun  = new G4GeneralParticleSource();
    
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


QuasiparticlePrimaryGeneratorAction::~QuasiparticlePrimaryGeneratorAction() {
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

 
void QuasiparticlePrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {

  /*
  if (fParticleGun->GetParticleDefinition() == G4Geantino::Definition()) {
    G4double selector = G4UniformRand();
    if (selector<0.53539) {
      fParticleGun->SetParticleDefinition(G4PhononTransSlow::Definition()); 
    } else if (selector<0.90217) {
      fParticleGun->SetParticleDefinition(G4PhononTransFast::Definition());
    } else {
      fParticleGun->SetParticleDefinition(G4PhononLong::Definition());
    }
  }
  */
  //  fParticleGun->SetParticleMomentumDirection(G4RandomDirection());
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



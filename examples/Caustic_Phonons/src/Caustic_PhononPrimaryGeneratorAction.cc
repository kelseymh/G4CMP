/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/


//20240110 Israel Hernandez -- Illinois Institute of Technology
//We include the GeneralParticle Source to control the distribution of the Particles
//The initial population is only to generate the caustic patterns
#include "Caustic_PhononPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4Geantino.hh"
#include "G4ParticleGun.hh"
#include "G4RandomDirection.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononTransSlow.hh"
#include "G4PhononLong.hh"
#include "G4SystemOfUnits.hh"

using namespace std;

Caustic_PhononPrimaryGeneratorAction::Caustic_PhononPrimaryGeneratorAction() {

fParticleGun  = new G4GeneralParticleSource();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


Caustic_PhononPrimaryGeneratorAction::~Caustic_PhononPrimaryGeneratorAction() {
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void Caustic_PhononPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
  //Defining 50% of the initial  population for TransFast and TransSlow phonons, this is only to distinguish the Phonons Caustics
  //Note.- This is not the initial density of States of the Substrate
  G4double selector = G4UniformRand();
  if (selector<0.5) {
  fParticleGun->SetParticleDefinition(G4PhononTransFast::Definition());

}
  else{
  fParticleGun->SetParticleDefinition(G4PhononTransSlow::Definition());
  // fParticleGun->SetParticleDefinition(G4PhononLong::Definition()); If you are interested in Longitudinal phonons, but their direction does not change too much.
  //You only need to uncomment and comment on the other fParticleGun.
  }
     fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

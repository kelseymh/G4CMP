
#include "ChargePrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4RandomDirection.hh"

#include "G4PhononTransFast.hh"
#include "G4PhononTransSlow.hh"
#include "G4PhononLong.hh"

#include "G4Electron.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"


using namespace std;

ChargePrimaryGeneratorAction::ChargePrimaryGeneratorAction()
{
  G4cout<<"\n\nChargePrimaryGeneratorAction::Constructor: running"<<endl;
  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);
  
  // default particle kinematic
  particleGun->SetParticleDefinition(G4CMPDriftElectron::Definition());
  particleGun->SetParticleMomentumDirection(G4ThreeVector(1,0,0));
  particleGun->SetParticlePosition(G4ThreeVector(0.0,0.0,0.0));
  particleGun->SetParticleEnergy(DBL_MIN);
}

ChargePrimaryGeneratorAction::~ChargePrimaryGeneratorAction()
{
  delete particleGun;
}

void ChargePrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  for(int i=0;i<10;i++)
  { 
    particleGun->SetParticlePosition(G4ThreeVector(0.0,0.0,0.0));    
    particleGun->SetParticleEnergy(1e-13*eV);  
    particleGun->SetParticleDefinition(G4CMPDriftElectron::Definition());
    particleGun->GeneratePrimaryVertex(anEvent);
    particleGun->SetParticleDefinition(G4CMPDriftHole::Definition());
    particleGun->GeneratePrimaryVertex(anEvent);
  }
}


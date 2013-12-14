
#include "Tst1PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4RandomDirection.hh"

#include "TPhononFast.hh"
#include "TPhononSlow.hh"
#include "LPhonon.hh"

#include "G4Electron.hh"
#include "DriftingElectron.hh"
#include "DriftingHole.hh"

#include "AlminumElectrodeSensitivity.hh"
#include "AlminumElectrodeHit.hh"

using namespace std;

Tst1PrimaryGeneratorAction::Tst1PrimaryGeneratorAction()
{
  G4cout<<"\n\nTst1PrimaryGeneratorAction::Constructor: running"<<endl;
  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);
  
 

  // default particle kinematic
  particleGun->SetParticleDefinition(DriftingElectron::Definition());
  particleGun->SetParticleMomentumDirection(G4ThreeVector(1,0,0));
  particleGun->SetParticlePosition(G4ThreeVector(0.0,0.0,0.0));
  particleGun->SetParticleEnergy(DBL_MIN);
  //  particleGun->SetParticlePolarization(*(new G4ThreeVector(0,1,1)));
}

Tst1PrimaryGeneratorAction::~Tst1PrimaryGeneratorAction()
{
  delete particleGun;
}

void Tst1PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  for(int i=0;i<10;i++)
  { 
    //particleGun->SetParticleDefinition(DriftingHole::Definition()); 
    //particleGun->SetParticleDefinition(LPhonon::Definition());
    
    particleGun->SetParticlePosition(G4ThreeVector(0.0,0.0,1.26*cm));    
    
    //particleGun->SetParticlePosition(G4ThreeVector(0.0,0.0,0.5*cm));    
    

    //particleGun->SetParticlePosition(G4ThreeVector(0.0,0.0,2.52*cm));    
     
//particleGun->SetParticlePosition(G4ThreeVector(0.0,0.0,-1.26*cm));    
    //particleGun->SetParticleEnergy(9.67e-6*eV);  
    

    particleGun->SetParticleEnergy(1e-13*eV);  
    //particleGun->SetParticleEnergy(1e-15*eV);

    //particleGun->SetParticleEnergy(1*eV);  
    

    //particleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1));
    

    //particleGun->GeneratePrimaryVertex(anEvent);

    //particleGun->SetParticleDefinition(DriftingHole::Definition());   
    particleGun->SetParticleDefinition(DriftingElectron::Definition());
    particleGun->GeneratePrimaryVertex(anEvent);
    }
  /*
  for(int i=0;i<10;i++)
     { 


       G4double selector = G4UniformRand();
       if(selector<0.53539) {
	 //ST++;
	 particleGun->SetParticleDefinition(TPhononSlow::PhononDefinition()); 
	 //particleGun->SetParticleMomentumDirection(LatticeManager2::mapKtoVDir(NULL,1,G4RandomDirection()));
       }
       else if(selector<0.90217) {
	 //FT++;
	 particleGun->SetParticleDefinition(TPhononFast::PhononDefinition());
	 //particleGun->SetParticleMomentumDirection(LatticeManager2::mapKtoVDir(NULL,2,G4RandomDirection()));
       }
       else {
	 // L++;
	 particleGun->SetParticleDefinition(TPhononFast::PhononDefinition());
	 //particleGun->SetParticleMomentumDirection(LatticeManager2::mapKtoVDir(NULL,0,G4RandomDirection()));
       }
       particleGun->SetParticleMomentumDirection(G4RandomDirection());
       particleGun->SetParticleEnergy(0.025*eV);
       particleGun->GeneratePrimaryVertex(anEvent);

  }
  */
 
  //G4cout<<"\n Fraction of L-Phonons: "<<L/1000<<"\n ST-Phonons: "<<ST/1000<<"\n FT-Phonons: "<<FT/1000<<"\n";
}


/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id: 028bb9edc3a4a6377c793f64eb4399f7ca05cfc7 $

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ChannelingPrimaryGeneratorAction.hh"

#include "ChannelingDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

#include "ChannelingPrimaryGeneratorMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ChannelingPrimaryGeneratorAction::ChannelingPrimaryGeneratorAction()
{
    fMessenger = new ChannelingPrimaryGeneratorMessenger(this);
    
    fParticleGun = new G4ParticleGun(1);
    
    fParticleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("proton"));
    
    fParticleGun->SetParticleEnergy(400. * CLHEP::GeV);
    
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,1.,0.));
    
    fParticleGun->SetParticlePosition(G4ThreeVector(0.,-1400. * CLHEP::centimeter,0.));
    
    fDivDistribution = "";
    
    fCutX = 0.;
    
    fCutY = 0.;
    
    fDivX = 0.;
    
    fDivY = 0.;
    
    fWidthX = 0.;
    
    fWidthY = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ChannelingPrimaryGeneratorAction::~ChannelingPrimaryGeneratorAction()
{
    delete fParticleGun;
    delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChannelingPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    //----------------------------------------
    // Function called at the beginning of an event
    //--------------------------------------
    
    G4ThreeVector vParticleMomentumDirection = G4ThreeVector(0.,1.,0.);
    
    G4double vRotationX = 0.;
    G4double vRotationY = 0.;
    
    G4bool bCutX = true;
    G4bool bCutY = true;
    
    do{
        if(fDivDistribution == "gauss"){
            vRotationX = CLHEP::RandGauss::shoot(0.,fDivX);
            vRotationY = CLHEP::RandGauss::shoot(0.,fDivY);
        }
        else{
            vRotationX = (G4UniformRand() - 0.5 ) * 2. * fDivX;
            vRotationY = (G4UniformRand() - 0.5 ) * 2. * fDivY;
        }
        if(fabs(vRotationX)<fCutX ||
           fCutX == 0.){
            bCutX = true;
        }
        else{
            bCutX = false;
        }
        if(fabs(vRotationY)<fCutY ||
           fCutY == 0.){
            bCutY = true;
        }
        else{
            bCutY = false;
        }
    }
    while(!(bCutX && bCutY));
    
    vParticleMomentumDirection = G4ThreeVector(0.,1.,0.).rotate(G4ThreeVector(0,0,1),vRotationX).unit();
    
    vParticleMomentumDirection = vParticleMomentumDirection.rotate(G4ThreeVector(0,1,0),vRotationY).unit();
    
    fParticleGun->SetParticleMomentumDirection(vParticleMomentumDirection);
    
    G4ThreeVector vPosition = fParticleGun->GetParticlePosition();
    
    vPosition.setX(CLHEP::RandGauss::shoot(0.,fWidthX));
    
    vPosition.setZ(CLHEP::RandGauss::shoot(0.,fWidthY));
    
    fParticleGun->SetParticlePosition(vPosition);
        
    fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


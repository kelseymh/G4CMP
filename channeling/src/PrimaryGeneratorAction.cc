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
//
// $Id: PrimaryGeneratorAction.cc,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

#include "PrimaryGeneratorActionMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
    fMessenger = new PrimaryGeneratorActionMessenger(this);
    
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

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fParticleGun;
    delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
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


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


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
    G4int vNumberOfParticles = 1;
    fParticleGun  = new G4ParticleGun(vNumberOfParticles);
    
    G4ParticleDefinition* vParticle = G4ParticleTable::GetParticleTable()->FindParticle("proton");
    
    fParticleGun->SetParticleDefinition(vParticle);
    
    G4ThreeVector vParticlePosition = G4ThreeVector(0.,-60.*cm,0.);
    
    fParticleGun->SetParticleEnergy(400.*GeV);
    
    fParticleGun->SetParticlePosition(vParticlePosition);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    //----------------------------------------
    // Function called at the beginning of an event
    //--------------------------------------
    
    G4ThreeVector vParticleMomentumDirection = G4ThreeVector(0.,1.,0.);

    G4bool bCollimated = false;

    if(!bCollimated){
        G4double vBeamDivergence = 15.0E-6 * radian;
        G4double vRotation = (G4UniformRand() - 0.5 ) * 2. * vBeamDivergence;
        vParticleMomentumDirection = G4ThreeVector(0.,1.,0.).rotate(G4ThreeVector(0,0,1),vRotation).unit();
    }
    
    fParticleGun->SetParticleMomentumDirection(vParticleMomentumDirection);
    
    fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


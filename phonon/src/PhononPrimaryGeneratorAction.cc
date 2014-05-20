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
/// \file exoticphysics/phonon/src/PhononPrimaryGeneratorAction.cc
/// \brief Implementation of the PhononPrimaryGeneratorAction class
//
// $Id$
//
// 20140519  Allow the user to specify phonon type by name in macro; if
//	     "geantino" is set, use random generator to select.

#include "PhononPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4Geantino.hh"
#include "G4ParticleGun.hh"
#include "G4RandomDirection.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononTransSlow.hh"
#include "G4PhononLong.hh"
#include "G4SystemOfUnits.hh"

using namespace std;

PhononPrimaryGeneratorAction::PhononPrimaryGeneratorAction() { 
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);   

  // default particle kinematics ("geantino" triggers random phonon choice)
  fParticleGun->SetParticleDefinition(G4Geantino::Definition());
  fParticleGun->SetParticleMomentumDirection(G4RandomDirection());
  fParticleGun->SetParticlePosition(G4ThreeVector(0.0,0.0,0.0));
  fParticleGun->SetParticleEnergy(0.0075*eV);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


PhononPrimaryGeneratorAction::~PhononPrimaryGeneratorAction() {
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

 
void PhononPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
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

  fParticleGun->SetParticleMomentumDirection(G4RandomDirection());
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



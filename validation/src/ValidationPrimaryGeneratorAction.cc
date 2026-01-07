/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file exoticphysics/validation/src/ValidationPrimaryGeneratorAction.cc
/// \brief Implementation of the ValidationPrimaryGeneratorAction class

#include "ValidationPrimaryGeneratorAction.hh"
#include "G4Event.hh"
#include "G4Geantino.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleGun.hh"
#include "G4RandomDirection.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononTransSlow.hh"
#include "G4PhononLong.hh"
#include "G4SystemOfUnits.hh"

using namespace std;

ValidationPrimaryGeneratorAction::ValidationPrimaryGeneratorAction() { 

  fParticleGun  = new G4GeneralParticleSource();
}


ValidationPrimaryGeneratorAction::~ValidationPrimaryGeneratorAction() {
  delete fParticleGun;
}

 
void ValidationPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {

  fParticleGun->GeneratePrimaryVertex(anEvent);
}

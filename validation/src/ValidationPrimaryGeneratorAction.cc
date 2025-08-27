/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file exoticphysics/validation/src/ValidationPrimaryGeneratorAction.cc
/// \brief Implementation of the ValidationPrimaryGeneratorAction class
//
// $Id: e75f788b103aef810361fad30f75077829192c13 $
//
// 20140519  Allow the user to specify phonon type by name in macro; if
//	     "geantino" is set, use random generator to select.

#include "ValidationPrimaryGeneratorAction.hh"

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

ValidationPrimaryGeneratorAction::ValidationPrimaryGeneratorAction() { 

  fParticleGun  = new G4GeneralParticleSource();
}


ValidationPrimaryGeneratorAction::~ValidationPrimaryGeneratorAction() {
  delete fParticleGun;
}

 
void ValidationPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {

  fParticleGun->GeneratePrimaryVertex(anEvent);
}



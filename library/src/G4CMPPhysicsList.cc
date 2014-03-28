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
/// \file exoticphysics/phonon/src/G4CMPPhysicsList.cc
/// \brief Implementation of the G4CMPPhysicsList class
//
// $Id$
//
// 20140328  Change to ModularPhysicsList, use new G4CMPPhysics for processes

#include "G4CMPPhysicsList.hh"
#include "G4CMPPhysics.hh"


// Constructor and destructor

G4CMPPhysicsList::G4CMPPhysicsList(G4int verbose) : G4VModularPhysicsList() {
  SetVerboseLevel(verbose);
  if (verbose) G4cout << "G4CMPPhysicsList::constructor" << G4endl;

  defaultCutValue = DBL_MIN;	// 100*mm;

  RegisterPhysics(new G4CMPPhysics);		// Phonons and charge-carriers
}

G4CMPPhysicsList::~G4CMPPhysicsList() {;}

// These values are used as the default production thresholds
// for the world volume.

void G4CMPPhysicsList::SetCuts() {
  SetCutsWithDefault();
}



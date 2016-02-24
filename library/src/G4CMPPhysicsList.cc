/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

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



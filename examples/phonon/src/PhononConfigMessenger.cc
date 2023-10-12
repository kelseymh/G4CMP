/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
// File:  PhononConfigMessenger.cc
//
// Description:	Macro command defitions to set user configuration in
//		PhononConfigManager.
//
// 20170816  Michael Kelsey

#include "PhononConfigMessenger.hh"
#include "PhononConfigManager.hh"
#include "G4UIcmdWithAString.hh"


// Constructor and destructor

PhononConfigMessenger::PhononConfigMessenger(PhononConfigManager* mgr)
  : G4UImessenger("/g4cmp/", "User configuration for G4CMP phonon example"),
    theManager(mgr), hitsCmd(0) {
  hitsCmd = CreateCommand<G4UIcmdWithAString>("HitsFile",
			      "Set filename for output of phonon hit locations");
}


PhononConfigMessenger::~PhononConfigMessenger() {
  delete hitsCmd; hitsCmd=0;
}


// Parse user input and add to configuration

void PhononConfigMessenger::SetNewValue(G4UIcommand* cmd, G4String value) {
  if (cmd == hitsCmd) theManager->SetHitOutput(value);
}

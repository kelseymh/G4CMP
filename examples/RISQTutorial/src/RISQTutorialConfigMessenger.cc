/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
// File:  RISQTutorialConfigMessenger.cc
//
// Description:	Macro command defitions to set user configuration in
//		RISQTutorialConfigManager.
//
// 20170816  Michael Kelsey

#include "RISQTutorialConfigMessenger.hh"
#include "RISQTutorialConfigManager.hh"
#include "G4UIcmdWithAString.hh"


// Constructor and destructor

RISQTutorialConfigMessenger::RISQTutorialConfigMessenger(RISQTutorialConfigManager* mgr)
  : G4UImessenger("/g4cmp/", "User configuration for G4CMP phonon example"),
    theManager(mgr), hitsCmd(0) {
  hitsCmd = CreateCommand<G4UIcmdWithAString>("HitsFile",
			      "Set filename for output of phonon hit locations");
}


RISQTutorialConfigMessenger::~RISQTutorialConfigMessenger() {
  delete hitsCmd; hitsCmd=0;
}


// Parse user input and add to configuration

void RISQTutorialConfigMessenger::SetNewValue(G4UIcommand* cmd, G4String value) {
  if (cmd == hitsCmd) theManager->SetHitOutput(value);
}

/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/


// 20241024 Israel Hernandez -- IIT, QSC and Fermilab

#include "Caustic_PhononConfigMessenger.hh"
#include "Caustic_PhononConfigManager.hh"
#include "G4UIcmdWithAString.hh"


// Constructor and destructor

Caustic_PhononConfigMessenger::Caustic_PhononConfigMessenger(Caustic_PhononConfigManager* mgr)
  : G4UImessenger("/g4cmp/", "User configuration for G4CMP phonon example"),
    theManager(mgr), hitsCmd(0) {
  hitsCmd = CreateCommand<G4UIcmdWithAString>("HitsFile",
			      "Set filename for output of phonon hit locations");

}


Caustic_PhononConfigMessenger::~Caustic_PhononConfigMessenger() {
  delete hitsCmd; hitsCmd=0;

}


// Parse user input and add to configuration

void Caustic_PhononConfigMessenger::SetNewValue(G4UIcommand* cmd, G4String value) {
  if (cmd == hitsCmd) theManager->SetHitOutput(value);

}

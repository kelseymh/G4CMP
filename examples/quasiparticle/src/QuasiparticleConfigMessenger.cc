/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
// File:  QuasiparticleConfigMessenger.cc
//
// Description:	Macro command defitions to set user configuration in
//		QuasiparticleConfigManager.
//
// 20170816  Michael Kelsey

#include "QuasiparticleConfigMessenger.hh"
#include "QuasiparticleConfigManager.hh"
#include "G4UIcmdWithAString.hh"


// Constructor and destructor

QuasiparticleConfigMessenger::
QuasiparticleConfigMessenger(QuasiparticleConfigManager* mgr)
  : G4UImessenger("/g4cmp/", "User configuration for G4CMP phonon example"),
    theManager(mgr), hitsCmd(0) {
  hitsCmd =
    CreateCommand<G4UIcmdWithAString>("HitsFile",
                                      "Set filename for output of phonon hit locations");
}


QuasiparticleConfigMessenger::~QuasiparticleConfigMessenger() {
  delete hitsCmd; hitsCmd=0;
}


// Parse user input and add to configuration

void QuasiparticleConfigMessenger::SetNewValue(G4UIcommand* cmd, G4String value) {
  if (cmd == hitsCmd) theManager->SetHitOutput(value);
}

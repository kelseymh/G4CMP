/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file ValidationConfigMessenger.cc
/// \brief	Macro command defitions to set user configuration in
///		ValidationConfigManager.

#include "ValidationConfigMessenger.hh"
#include "ValidationConfigManager.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"

// Constructor and destructor

ValidationConfigMessenger::ValidationConfigMessenger(ValidationConfigManager* mgr)
  : G4UImessenger("/g4cmp/", "User configuration for G4CMP phonon example"),
    theManager(mgr), hitsCmd(0) {
  hitsCmd = CreateCommand<G4UIcmdWithAString>("HitsFile",
                                              "Set filename for output of phonon hit locations");
  geometryCmd = CreateCommand<G4UIcmdWithAnInteger>("GeometryID",
                                                    "Set the validation geometry ID for use");
  stepFileCmd = CreateCommand<G4UIcmdWithAString>("StepFile",
                                                  "Set filename for step output file (validations).");
}


ValidationConfigMessenger::~ValidationConfigMessenger() {
  delete hitsCmd; hitsCmd=0;
  delete geometryCmd; geometryCmd=0;
  delete stepFileCmd; stepFileCmd=0;
}


// Parse user input and add to configuration

void ValidationConfigMessenger::SetNewValue(G4UIcommand* cmd, G4String value) {
  if (cmd == hitsCmd) theManager->SetHitOutput(value);
  if (cmd == geometryCmd) theManager->SetGeometryID(StoI(value));
  if (cmd == stepFileCmd) theManager->SetStepFile(value);
}

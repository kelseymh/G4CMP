/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
// File:  ChargeConfigMessenger.cc
//
// Description:	Macro command defitions to set user configuration in
//		ChargeConfigManager.
//
// 20170816  Michael Kelsey

#include "ChargeConfigMessenger.hh"
#include "ChargeConfigManager.hh"
#include "G4Tokenizer.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"


// Constructor and destructor

ChargeConfigMessenger::ChargeConfigMessenger(ChargeConfigManager* mgr)
  : G4UImessenger("/g4cmp/", "User configuration for G4CMP phonon example"),
    theManager(mgr), voltageCmd(0), escaleCmd(0), fileCmd(0), hitsCmd(0),
    millerCmd(0) {
  voltageCmd = CreateCommand<G4UIcmdWithADoubleAndUnit>("voltage",
				"Set voltage for uniform electric field");
  voltageCmd->SetUnitCategory("Electric potential");

  escaleCmd = CreateCommand<G4UIcmdWithADouble>("scaleEPot",
		"Set a scale factor for voltages in EPot electric field file");

  fileCmd = CreateCommand<G4UIcmdWithAString>("EPotFile",
			      "Set filename for non-uniform electric field");

  hitsCmd = CreateCommand<G4UIcmdWithAString>("HitsFile",
			      "Set filename for output of phonon hit locations");

  millerCmd = CreateCommand<G4UIcmdWithAString>("orientation",
	"Lattice orientation (Miller indices h, k, l in direct basis)");
  millerCmd->SetGuidance("Orientation aligns with Z axis of G4VSolid");
}


ChargeConfigMessenger::~ChargeConfigMessenger() {
  delete voltageCmd; voltageCmd=0;
  delete escaleCmd; escaleCmd=0;
  delete fileCmd; fileCmd=0;
  delete hitsCmd; hitsCmd=0;
  delete millerCmd; millerCmd=0;
}


// Parse user input and add to configuration

void ChargeConfigMessenger::SetNewValue(G4UIcommand* cmd, G4String value) {
  if (cmd == fileCmd) theManager->SetEPotFile(value);
  if (cmd == hitsCmd) theManager->SetHitOutput(value);

  if (cmd == voltageCmd)
    theManager->SetVoltage(voltageCmd->GetNewDoubleValue(value));

  if (cmd == escaleCmd)
    theManager->SetEPotScale(escaleCmd->GetNewDoubleValue(value));

  if (cmd == millerCmd) {		// Special, takes three integer args
    G4Tokenizer split(value);
    theManager->SetMillerOrientation(StoI(split()),StoI(split()),StoI(split()));
  }
}

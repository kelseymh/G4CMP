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
// 20260422  G4CMP-597 -- Make fileCmd optional and provide default value.
// 20260513  G4CMP-604 -- Add surface property UI commands to charge example.

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
    millerCmd(0), absorbCmd(0), specularCmd(0) {
  voltageCmd = CreateCommand<G4UIcmdWithADoubleAndUnit>("voltage",
				"Set voltage for uniform electric field");
  voltageCmd->SetUnitCategory("Electric potential");

  escaleCmd = CreateCommand<G4UIcmdWithADouble>("scaleEPot",
		"Set a scale factor for voltages in EPot electric field file");

  fileCmd = CreateCommand<G4UIcmdWithAString>("EPotFile",
			      "Set filename for non-uniform electric field");
  // Make the filename argument optional with empty string as default
  fileCmd->SetParameterName("filename", true);
  fileCmd->SetDefaultValue("");

  hitsCmd = CreateCommand<G4UIcmdWithAString>("HitsFile",
			      "Set filename for output of phonon hit locations");

  millerCmd = CreateCommand<G4UIcmdWithAString>("orientation",
	"Lattice orientation (Miller indices h, k, l in direct basis)");
  millerCmd->SetGuidance("Orientation aligns with Z axis of G4VSolid");

  absorbCmd = CreateCommand<G4UIcmdWithADouble>("chargeAbsorb",
		"Set probability (0-1) for charge carriers to be absorbed at detector surfaces (top, bottom, sidewalls)");

  specularCmd = CreateCommand<G4UIcmdWithADouble>("specularReflect",
		"Set probability (0-1) for specular reflection of charge carriers at detector surfaces (top, bottom, sidewalls)");
}


ChargeConfigMessenger::~ChargeConfigMessenger() {
  delete voltageCmd; voltageCmd=0;
  delete escaleCmd; escaleCmd=0;
  delete fileCmd; fileCmd=0;
  delete hitsCmd; hitsCmd=0;
  delete millerCmd; millerCmd=0;
  delete absorbCmd; absorbCmd=0;
  delete specularCmd; specularCmd=0;
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

  if (cmd == absorbCmd)
    theManager->SetChargeAbsorbProb(absorbCmd->GetNewDoubleValue(value));

  if (cmd == specularCmd)
    theManager->SetSpecularReflectProb(specularCmd->GetNewDoubleValue(value));
}

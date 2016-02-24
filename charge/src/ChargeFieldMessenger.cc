/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#include "ChargeFieldMessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "ChargeEMField.hh"
#include "G4UIcommandTree.hh"
#include "G4UImanager.hh"

ChargeFieldMessenger::ChargeFieldMessenger(ChargeEMField* field):
                      ChargeField(field),localFieldDir(false)
{
  // See if input path has already been registered
  G4UImanager* UIman = G4UImanager::GetUIpointer();
  G4UIcommand* foundPath = UIman->GetTree()->FindPath("/field/");
  if (foundPath) fieldDir = dynamic_cast<G4UIdirectory*>(foundPath);

  if (!fieldDir) {		// Create local deletable directory
    localFieldDir = true;
    fieldDir = new G4UIdirectory("/field/");
    fieldDir->SetGuidance("EM Field commands");
  }

  fieldTree = UIman->GetTree()->FindCommandTree("/field/");	// For printing

  buildcmd = new G4UIcmdWithoutParameter("/field/Initialize", this);
  buildcmd->SetGuidance("Construct electric field.");

  constFieldcmd = new G4UIcmdWithABool("/field/UseConstantField",this);
  constFieldcmd->SetGuidance("Use a constant electric field instead of interpolating.");
  constFieldcmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  constFieldMagcmd = new G4UIcmdWithADoubleAndUnit("/field/ConstantFieldMag", this);
  constFieldMagcmd->SetGuidance("Set strength of constant electric field. [158 volt/m]");
  constFieldMagcmd->SetParameterName("Electric Field Magntitude", false);
  constFieldMagcmd->SetUnitCandidates("volt/m volt/cm");
  constFieldMagcmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  constFieldDircmd = new G4UIcmdWith3Vector("/field/ConstantFieldDir", this);
  constFieldDircmd->SetGuidance("Set direction of constant electric field. [(0, 0, 1)]");
  constFieldDircmd->SetParameterName("Electric Field X Direction",
                                     "Electric Field Y Direction",
                                     "Electric Field Z Direction", false);
  constFieldDircmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fieldFileNamecmd = new G4UIcmdWithAString("/field/EPotFileName", this);
  fieldFileNamecmd->SetGuidance("EPot file name for field interpolation.");
  fieldFileNamecmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

ChargeFieldMessenger::~ChargeFieldMessenger()
{
  if (localFieldDir) delete fieldDir;
  if (localFieldDir) delete fieldTree;
  delete buildcmd;
  delete constFieldcmd;
  delete constFieldMagcmd;
  delete constFieldDircmd;
  delete fieldFileNamecmd;
}

void ChargeFieldMessenger::SetNewValue(G4UIcommand* command, G4String NewValue)
{
    if (command == buildcmd)
      ChargeField->Build();
    else if (command == fieldFileNamecmd)
      ChargeField->SetEPotFileName(NewValue);
    else if (command == constFieldcmd)
      ChargeField->UseConstantField(constFieldcmd->GetNewBoolValue(NewValue));
    else if (command == constFieldMagcmd)
      ChargeField->SetFieldMagnitude(constFieldMagcmd->GetNewDoubleValue(NewValue));
    else if (command == constFieldDircmd)
      ChargeField->SetFieldDirection(constFieldDircmd->GetNew3VectorValue(NewValue));
}

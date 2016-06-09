/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file hadronic/Hadr01/src/ChannelingStackingMessenger.cc
/// \brief Implementation of the ChannelingStackingMessenger class
//
// $Id: 6ce68fdf9bcb0a144a1d583b595642fc9be52f22 $
//

#include "ChannelingStackingMessenger.hh"
#include "ChannelingStackingAction.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ChannelingStackingMessenger::ChannelingStackingMessenger(ChannelingStackingAction* stack)
:fStackAction(stack)
{
  fKillCmd = new G4UIcmdWithABool("/mystack/KillAllSecondaries",this);
  fKillCmd->SetGuidance("  Choice : true false");
  fKillCmd->SetParameterName("choice",true);
  fKillCmd->SetDefaultValue(false);

  fKCmd = new G4UIcmdWithAString("/mystack/Kill", this);
  fKCmd->SetGuidance("Kill secondary particles of defined type");
  fKCmd->SetParameterName("ch", true);
  fKCmd->SetDefaultValue("none");
  fKCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ChannelingStackingMessenger::~ChannelingStackingMessenger()
{
  delete fKillCmd;
  delete fKCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChannelingStackingMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{     
  if(command == fKillCmd) {
    fStackAction->SetKillStatus(fKillCmd->GetNewBoolValue(newValue));
  } else if(command == fKCmd) {
    fStackAction->SetKill(newValue);               
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

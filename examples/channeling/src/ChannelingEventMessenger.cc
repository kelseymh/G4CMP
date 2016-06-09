/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file analysis/A01/src/ChannelingEventMessenger.cc
/// \brief Implementation of the ChannelingEventMessenger class
//
// $Id: 0bd6e89b9f3bc221e15166d7b82c65ed54e35d03 $
// --------------------------------------------------------------
//

#include "ChannelingEventMessenger.hh"
#include "ChannelingEventAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ios.hh"

ChannelingEventMessenger::ChannelingEventMessenger(ChannelingEventAction * mpga)
:fTarget (mpga)
{
    fMyDirectory = new G4UIdirectory("/myevt/");
    
    fFileNameCmd = new G4UIcmdWithAString("/myevt/filename",this);
    fFileNameCmd->SetGuidance("Filename for output.");
    fFileNameCmd->SetParameterName("evtfilename",true);
    fFileNameCmd->SetDefaultValue("noname");
    
    fVerboseCmd = new G4UIcmdWithAnInteger("/myevt/verbose",this);
    fVerboseCmd->SetGuidance("Verbose level for each event.");
    fVerboseCmd->SetGuidance(" Event summary will be displayed for every 'level' events.");
    fVerboseCmd->SetParameterName("level",true);
    fVerboseCmd->SetRange("level>=0");
    fVerboseCmd->SetDefaultValue(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ChannelingEventMessenger::~ChannelingEventMessenger()
{
    delete fVerboseCmd;
    delete fFileNameCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ChannelingEventMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
    if( command==fVerboseCmd ){
        fTarget->SetVerbose(fVerboseCmd->GetNewIntValue(newValue));
    }
    
    if( command==fFileNameCmd ){
        fTarget->SetFileName(newValue);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4String ChannelingEventMessenger::GetCurrentValue(G4UIcommand * command)
{
    G4String cv;
    if( command==fVerboseCmd ){
        cv = fVerboseCmd->ConvertToString(fTarget->GetVerbose());
    }
    
    if( command==fFileNameCmd ){
        cv = fTarget->GetFileName();
    }
    
    return cv;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

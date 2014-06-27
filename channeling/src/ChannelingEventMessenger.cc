//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file analysis/A01/src/ChannelingEventMessenger.cc
/// \brief Implementation of the ChannelingEventMessenger class
//
// $Id$
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

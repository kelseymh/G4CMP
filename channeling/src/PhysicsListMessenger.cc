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
/// \file analysis/A01/src/PhysicsListMessenger.cc
/// \brief Implementation of the PhysicsListMessenger class
//
// $Id$
// --------------------------------------------------------------
//

#include "PhysicsListMessenger.hh"
#include "PhysicsList.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ios.hh"

PhysicsListMessenger::PhysicsListMessenger(PhysicsList * mpga)
:fTarget (mpga)
{
    fMyDirectory = new G4UIdirectory("/myproc/");
    
    fFileNameCmd = new G4UIcmdWithAString("/myproc/filename",this);
    fFileNameCmd->SetGuidance("Filename for output.");
    fFileNameCmd->SetParameterName("evtfilename",true);
    fFileNameCmd->SetDefaultValue("noname");

    fScatteringType = new G4UIcmdWithAString("/myproc/scattering",this);
    fScatteringType->SetGuidance("Scattering type.");
    fScatteringType->SetParameterName("scattype",true);
    fScatteringType->SetDefaultValue("ss");

    fChannelingCmd = new G4UIcmdWithABool("/myproc/channeling",this);
    fChannelingCmd->SetGuidance("Enable channeling");
    fChannelingCmd->SetParameterName("channeling",true);
    fChannelingCmd->SetDefaultValue(1);
    
    fWrapperCmd = new G4UIcmdWithABool("/myproc/wrapper",this);
    fWrapperCmd->SetGuidance("Enable wrapper");
    fWrapperCmd->SetParameterName("wrapper",true);
    fWrapperCmd->SetDefaultValue(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PhysicsListMessenger::~PhysicsListMessenger()
{
    delete fFileNameCmd;
    delete fScatteringType;
    delete fChannelingCmd;
    delete fWrapperCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsListMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
    if( command==fFileNameCmd ){
        fTarget->SetFileName(newValue);
    }
    if( command==fScatteringType ){
        fTarget->SetScatteringType(newValue);
    }
    if( command==fChannelingCmd ){
        fTarget->EnableChanneling(fChannelingCmd->GetNewBoolValue(newValue));
    }
    if( command==fWrapperCmd ){
        fTarget->EnableWrapper(fWrapperCmd->GetNewBoolValue(newValue));
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4String PhysicsListMessenger::GetCurrentValue(G4UIcommand * command)
{
    G4String cv;
    if( command==fFileNameCmd ){
        cv = fTarget->GetFileName();
    }
    if( command==fScatteringType ){
        cv = fTarget->GetScatteringType();
    }
    if( command==fChannelingCmd ){
        cv = fChannelingCmd->ConvertToString(fTarget->GetChannelingState());
    }
    if( command==fChannelingCmd ){
        cv = fChannelingCmd->ConvertToString(fTarget->GetWrapperState());
    }
    
    return cv;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

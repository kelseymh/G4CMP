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
/// \file analysis/A01/src/PrimaryGeneratorActionMessenger.cc
/// \brief Implementation of the PrimaryGeneratorActionMessenger class
//
// $Id$
// --------------------------------------------------------------
//

#include "PrimaryGeneratorActionMessenger.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4ios.hh"

PrimaryGeneratorActionMessenger::PrimaryGeneratorActionMessenger(PrimaryGeneratorAction * mpga)
:fTarget (mpga)
{
    fDivergenceX = new G4UIcmdWithADoubleAndUnit("/gun/setDivX",this);
    fDivergenceX->SetGuidance("Set beam divergence Y.");
    fDivergenceX->SetParameterName("beamdivx",true);
    fDivergenceX->SetDefaultValue(0.);
    fDivergenceX->SetDefaultUnit("rad");

    fDivergenceY = new G4UIcmdWithADoubleAndUnit("/gun/setDivY",this);
    fDivergenceY->SetGuidance("Set beam divergence X.");
    fDivergenceY->SetParameterName("beamdivy",true);
    fDivergenceY->SetDefaultValue(0.);
    fDivergenceY->SetDefaultUnit("rad");

    fWidthX = new G4UIcmdWithADoubleAndUnit("/gun/setWidthX",this);
    fWidthX->SetGuidance("Set beam width X.");
    fWidthX->SetParameterName("beamwidthx",true);
    fWidthX->SetDefaultValue(0.);
    fWidthX->SetDefaultUnit("mm");

    fWidthY = new G4UIcmdWithADoubleAndUnit("/gun/setWidthY",this);
    fWidthY->SetGuidance("Set beam width Y.");
    fWidthY->SetParameterName("beamwidthy",true);
    fWidthY->SetDefaultValue(0.);
    fWidthY->SetDefaultUnit("mm");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PrimaryGeneratorActionMessenger::~PrimaryGeneratorActionMessenger()
{
    delete fDivergenceX;
    delete fDivergenceY;
    delete fWidthX;
    delete fWidthY;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PrimaryGeneratorActionMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
    if(command==fDivergenceX ){
        fTarget->SetBeamDivergenceX(fDivergenceX->GetNewDoubleValue(newValue));
    }
    if(command==fDivergenceY ){
        fTarget->SetBeamDivergenceY(fDivergenceY->GetNewDoubleValue(newValue));
    }
    if(command==fWidthX ){
        fTarget->SetBeamWidthX(fWidthX->GetNewDoubleValue(newValue));
    }
    if(command==fWidthY ){
        fTarget->SetBeamWidthY(fWidthY->GetNewDoubleValue(newValue));
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4String PrimaryGeneratorActionMessenger::GetCurrentValue(G4UIcommand * command)
{
    G4String cv;
    if( command==fDivergenceX ){
        cv = fDivergenceX->ConvertToString(fTarget->GetBeamDivergenceX(),"rad");
    }
    if( command==fDivergenceY ){
        cv = fDivergenceY->ConvertToString(fTarget->GetBeamDivergenceY(),"rad");
    }
    if( command==fWidthX ){
        cv = fWidthX->ConvertToString(fTarget->GetBeamWidthX(),"rad");
    }
    if( command==fWidthY ){
        cv = fWidthY->ConvertToString(fTarget->GetBeamWidthY(),"rad");
    }
    
    return cv;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

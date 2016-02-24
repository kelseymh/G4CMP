/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file electromagnetic/TestEm5/src/ChannelingStepMessenger.cc
/// \brief Implementation of the ChannelingStepMessenger class
//
// $Id: 42e174909e17dbc65e38a3ab6385125b167fadfb $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ChannelingStepMessenger.hh"

#include "ChannelingStepLimiter.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ChannelingStepMessenger::ChannelingStepMessenger(ChannelingStepLimiter* stepM)
:fChannelingStepLimiter(stepM)
{ 
  fChannelingStepLimiterCmd = new G4UIcmdWithADoubleAndUnit("/testem/stepMax",this);
  fChannelingStepLimiterCmd->SetGuidance("Set max allowed step length");
  fChannelingStepLimiterCmd->SetParameterName("mxStep",false);
  fChannelingStepLimiterCmd->SetRange("mxStep>0.");
  fChannelingStepLimiterCmd->SetUnitCategory("Length");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ChannelingStepMessenger::~ChannelingStepMessenger()
{
  delete fChannelingStepLimiterCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChannelingStepMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if (command == fChannelingStepLimiterCmd)
    { fChannelingStepLimiter->SetMaxStep(fChannelingStepLimiterCmd->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

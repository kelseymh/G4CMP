#include "FETMessenger.hh"
#include "G4RunManager.hh"
#include "FET.hh"
#include "G4UIdirectory.hh"

FETMessenger::FETMessenger(FET* fet):target(fet)
{    
    fetDir = new G4UIdirectory("/FET/");
    fetDir->SetGuidance("FETSim commands");

    cmd = new G4UIcmdWithoutParameter("/FET/FETSim",this);
    cmd->SetGuidance("Run FETSim on most recent run");
}
FETMessenger::~FETMessenger()
{
    delete cmd;
    delete fetDir;
}

void FETMessenger::SetNewValue(G4UIcommand* command, G4String NewValue)
{
    if(command == cmd)
        target->Run();
}

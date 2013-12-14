#include "FETMessenger1.hh"
#include "G4RunManager.hh"
#include "FET1.hh"
#include "G4UIdirectory.hh"

FETMessenger1::FETMessenger1(FET1* fet):target(fet)
{    
    fetDir = new G4UIdirectory("/FET/");
    fetDir->SetGuidance("FETSim commands");

    cmd = new G4UIcmdWithoutParameter("/FET/FETSim",this);
    cmd->SetGuidance("Run FETSim on most recent run");
}
FETMessenger1::~FETMessenger1()
{
    delete cmd;
    delete fetDir;
}

void FETMessenger1::SetNewValue(G4UIcommand* command, G4String NewValue)
{
    if(command == cmd)
        target->Run();
}

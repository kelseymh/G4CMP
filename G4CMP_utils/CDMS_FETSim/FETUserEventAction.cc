#include"FETUserEventAction.hh" 
#include"G4EventManager.hh"

void FETUserEventAction::EndOfEventAction(const G4Event*)
{
    fpEventManager->KeepTheCurrentEvent();
}

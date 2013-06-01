#include "FETUserRunAction.hh"

FETUserRunAction::FETUserRunAction()
{;}

FETUserRunAction::~FETUserRunAction()
{;}
void FETUserRunAction::EndOfRunAction(const G4Run* aRun)
{
    G4cout << aRun->GetEventVector()->size() << G4endl;
    FET* FetSim = new FET(aRun);
    delete FetSim;
}

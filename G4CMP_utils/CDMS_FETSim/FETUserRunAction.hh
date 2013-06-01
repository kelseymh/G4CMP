#ifndef FET_User_Run_Action_hh
#define FET_User_Run_Action_hh
#include "FET.hh"
#include "G4UserRunAction.hh"
#include "G4Run.hh"

class FETUserRunAction : public G4UserRunAction
{
    public:
        FETUserRunAction();
        ~FETUserRunAction();

        void BeginOfRunAction(const G4Run* aRun) {}
        void EndOfRunAction(const G4Run* aRun);
};

#endif

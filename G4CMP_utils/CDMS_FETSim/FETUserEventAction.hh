#ifndef FET_USER_ACTION_hh
#define FET_USER_ACTION_hh

#include"globals.hh" 
#include "G4UserEventAction.hh"

class G4EventManager;
class FETUserEventAction : public G4UserEventAction
{
    public:
        FETUserEventAction() {;}
        ~FETUserEventAction() {;}
        void BeginofEventAction(const G4Event*) {;}
        void EndOfEventAction(const G4Event*);
};
#endif

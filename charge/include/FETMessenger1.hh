#ifndef FETMESSENGER_HH 
#define FETMESSENGER_HH 1

class FET1;
class G4UIdirectory;

#include "G4UImessenger.hh"
#include "globals.hh"
#include "G4UIcmdWithoutParameter.hh"

class FETMessenger1: public G4UImessenger
{
  public:
    FETMessenger1(FET1* fet);
    ~FETMessenger1();
  public:
    void SetNewValue(G4UIcommand* command, G4String NewValue);
  private:
    FET1* target;
  private:
    G4UIdirectory* fetDir;
    G4UIcmdWithoutParameter* cmd;
};

#endif

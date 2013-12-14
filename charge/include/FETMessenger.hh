#ifndef FETMESSENGER_HH 
#define FETMESSENGER_HH 1

class FET;
class G4UIdirectory;

#include "G4UImessenger.hh"
#include "globals.hh"
#include "G4UIcmdWithoutParameter.hh"

class FETMessenger: public G4UImessenger
{
  public:
    FETMessenger(FET* fet);
    ~FETMessenger();
  public:
    void SetNewValue(G4UIcommand* command, G4String NewValue);
  private:
    FET* target;
  private:
    G4UIdirectory* fetDir;
    G4UIcmdWithoutParameter* cmd;
};

#endif

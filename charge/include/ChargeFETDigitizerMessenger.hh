#ifndef CHARGEFETMESSENGER_HH
#define CHARGEFETMESSENGER_HH 1

class ChargeFETDigitizerModule;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;

#include "G4UImessenger.hh"
#include "globals.hh"

class ChargeFETDigitizerMessenger : public G4UImessenger
{
  public:
    ChargeFETDigitizerMessenger(ChargeFETDigitizerModule* digitizer);
    ~ChargeFETDigitizerMessenger();
  public:
    void SetNewValue(G4UIcommand* command, G4String NewValue);
  private:
    ChargeFETDigitizerModule* fet;
  private:
    G4UIdirectory*             fetDir;
    G4UIcmdWithAString*        SetOutputFileCmd;
    G4UIcmdWithoutParameter*   GetOutputFileCmd;
    G4UIcmdWithAString*        SetTemplateFileCmd;
    G4UIcmdWithoutParameter*   GetTemplateFileCmd;
    G4UIcmdWithAString*        SetRamoFileDirCmd;
    G4UIcmdWithoutParameter*   GetRamoFileCmd;
    G4UIcmdWithAnInteger*      SetNumChanCmd;
    G4UIcmdWithoutParameter*   GetNumChanCmd;
    G4UIcmdWithAnInteger*      SetTimeBinCmd;
    G4UIcmdWithoutParameter*   GetTimeBinCmd;
    G4UIcmdWithADoubleAndUnit* SetDecayTimeCmd;
    G4UIcmdWithoutParameter*   GetDecayTimeCmd;
    G4UIcmdWithADoubleAndUnit* SetUnitTimeCmd;
    G4UIcmdWithoutParameter*   GetUnitTimeCmd;
    G4UIcmdWithADoubleAndUnit* SetPreTrigCmd;
    G4UIcmdWithoutParameter*   GetPreTrigCmd;
    G4UIcmdWithADoubleAndUnit* SetTemplateEnergyCmd;
    G4UIcmdWithoutParameter*   GetTemplateEnergyCmd;
};

#endif

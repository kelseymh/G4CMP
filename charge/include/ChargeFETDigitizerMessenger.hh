/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef CHARGEFETMESSENGER_HH
#define CHARGEFETMESSENGER_HH 1

#include "G4UImessenger.hh"

class ChargeFETDigitizerModule;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithABool;

class ChargeFETDigitizerMessenger : public G4UImessenger
{
  public:
    ChargeFETDigitizerMessenger(ChargeFETDigitizerModule* digitizer);
    ~ChargeFETDigitizerMessenger();
    void SetNewValue(G4UIcommand* command, G4String NewValue);
  private:
    ChargeFETDigitizerModule*  fet;
    G4UIdirectory*             fetDir;
    G4UIcmdWithABool*          EnableFETSimCmd;
    G4UIcmdWithABool*          DisableFETSimCmd;
    G4UIcmdWithoutParameter*   GetEnabledStateCmd;
    G4UIcmdWithAString*        SetOutputFileCmd;
    G4UIcmdWithoutParameter*   GetOutputFileCmd;
    G4UIcmdWithAString*        SetConfigFileCmd;
    G4UIcmdWithoutParameter*   GetConfigFileCmd;
    G4UIcmdWithAString*        SetTemplateFileCmd;
    G4UIcmdWithoutParameter*   GetTemplateFileCmd;
    G4UIcmdWithAString*        SetRamoFileDirCmd;
    G4UIcmdWithoutParameter*   GetRamoFileDirCmd;
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
    G4UIcmdWithoutParameter*   UpdateCmd;
};

#endif

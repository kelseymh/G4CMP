/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#include "ChargeFETDigitizerMessenger.hh"
#include "ChargeFETDigitizerModule.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"

ChargeFETDigitizerMessenger::ChargeFETDigitizerMessenger(
                          ChargeFETDigitizerModule* digitizer) : fet(digitizer)
{    
  fetDir = new G4UIdirectory("/g4cmp/FETSim/");
  fetDir->SetGuidance("FETSim commands");

  EnableFETSimCmd = new G4UIcmdWithABool("/g4cmp/FETSim/EnableFETSim",this);
  EnableFETSimCmd->SetGuidance("Enable FET simulation during run.");

  DisableFETSimCmd = new G4UIcmdWithABool("/g4cmp/FETSim/DisableFETSim",this);
  DisableFETSimCmd->SetGuidance("Disable FET simulation during run.");

  GetEnabledStateCmd = new G4UIcmdWithoutParameter("/g4cmp/FETSim/GetEnabledState",this);
  GetEnabledStateCmd->SetGuidance("Report whether FET simulation is enabled.");

  SetOutputFileCmd = new G4UIcmdWithAString("/g4cmp/FETSim/SetOutputFile",this);
  SetOutputFileCmd->SetGuidance("Set path to FET output file.");

  GetOutputFileCmd = new G4UIcmdWithoutParameter("/g4cmp/FETSim/GetOutputFile",this);
  GetOutputFileCmd->SetGuidance("Current path to FET output file.");

  SetConfigFileCmd = new G4UIcmdWithAString("/g4cmp/FETSim/SetConfigFile",this);
  SetConfigFileCmd->SetGuidance("Set path to FET config file.");

  GetConfigFileCmd = new G4UIcmdWithoutParameter("/g4cmp/FETSim/GetConfigFile",this);
  GetConfigFileCmd->SetGuidance("Current path to FET config file.");

  SetTemplateFileCmd = new G4UIcmdWithAString("/g4cmp/FETSim/SetTemplateFile",this);
  SetTemplateFileCmd->SetGuidance("Set path to FET template file.");

  GetTemplateFileCmd = new G4UIcmdWithoutParameter("/g4cmp/FETSim/GetTemplateFile",this);
  GetTemplateFileCmd->SetGuidance("Current path to FET template file.");

  SetRamoFileDirCmd = new G4UIcmdWithAString("/g4cmp/FETSim/SetRamoFileDir",this);
  SetRamoFileDirCmd->SetGuidance("Set path to Ramo potential files.");

  GetRamoFileDirCmd = new G4UIcmdWithoutParameter("/g4cmp/FETSim/GetRamoFileDir",this);
  GetRamoFileDirCmd->SetGuidance("Current path to Ramo potential files.");

  SetNumChanCmd = new G4UIcmdWithAnInteger("/g4cmp/FETSim/SetNumberOfChannels",this);
  SetNumChanCmd->SetGuidance("Number of FET channels in detector needs to match number of Ramo potential files.");

  GetNumChanCmd = new G4UIcmdWithoutParameter("/g4cmp/FETSim/GetNumberOfChannels",this);
  GetNumChanCmd->SetGuidance("Number of FET channels");

  SetTimeBinCmd = new G4UIcmdWithAnInteger("/g4cmp/FETSim/SetNumberOfBins",this);
  SetTimeBinCmd->SetGuidance("Set number of digitizer bins");

  GetTimeBinCmd = new G4UIcmdWithoutParameter("/g4cmp/FETSim/GetNumberOfBins",this);
  GetTimeBinCmd->SetGuidance("Number of digitizer bins");

  SetDecayTimeCmd = new G4UIcmdWithADoubleAndUnit("/g4cmp/FETSim/SetDecayTime",this);
  SetDecayTimeCmd->SetGuidance("Pulse decay time (if not using templates)");

  GetDecayTimeCmd = new G4UIcmdWithoutParameter("/g4cmp/FETSim/GetDecayTime",this);
  GetDecayTimeCmd->SetGuidance("Get pulse decay time (if not using templates)");

  SetUnitTimeCmd = new G4UIcmdWithADoubleAndUnit("/g4cmp/FETSim/SetUnitTime",this);
  SetUnitTimeCmd->SetGuidance("Set dt for pulse bins (if not using templates)");

  GetUnitTimeCmd = new G4UIcmdWithoutParameter("/g4cmp/FETSim/GetUnitTime",this);
  GetUnitTimeCmd->SetGuidance("Get dt for pulse bins (if not using templates)");

  SetPreTrigCmd = new G4UIcmdWithADoubleAndUnit("/g4cmp/FETSim/SetPreTriggerTime",this);
  SetPreTrigCmd->SetGuidance("Set pre-trigger time for pulse (if not using templates)");

  GetPreTrigCmd = new G4UIcmdWithoutParameter("/g4cmp/FETSim/GetPreTriggerTime",this);
  GetPreTrigCmd->SetGuidance("Get pre-trigger time for pulse (if not using templates)");

  UpdateCmd = new G4UIcmdWithoutParameter("/g4cmp/FETSim/Update",this);
  UpdateCmd->SetGuidance("Must manually udpate FETSim after changing parameters.");
}

ChargeFETDigitizerMessenger::~ChargeFETDigitizerMessenger()
{
    delete fetDir;
    delete EnableFETSimCmd;
    delete DisableFETSimCmd;
    delete GetEnabledStateCmd;
    delete SetOutputFileCmd;
    delete GetOutputFileCmd;
    delete SetConfigFileCmd;
    delete GetConfigFileCmd;
    delete SetTemplateFileCmd;
    delete GetTemplateFileCmd;
    delete SetRamoFileDirCmd;
    delete GetRamoFileDirCmd;
    delete SetNumChanCmd;
    delete GetNumChanCmd;
    delete SetTimeBinCmd;
    delete GetTimeBinCmd;
    delete SetDecayTimeCmd;
    delete GetDecayTimeCmd;
    delete SetUnitTimeCmd;
    delete GetUnitTimeCmd;
    delete SetPreTrigCmd;
    delete GetPreTrigCmd;
    delete UpdateCmd;
}

void ChargeFETDigitizerMessenger::SetNewValue(G4UIcommand* command, G4String NewValue)
{
  if (command == EnableFETSimCmd)
    fet->EnableFETSim();
  else if (command == DisableFETSimCmd)
    fet->DisableFETSim();
  else if (command == SetOutputFileCmd)
    fet->SetOutputFile(NewValue);
  else if (command == GetOutputFileCmd)
    fet->GetOutputFile();
  else if (command == SetConfigFileCmd)
    fet->SetConfigFilename(NewValue);
  else if (command == GetConfigFileCmd)
    fet->GetConfigFilename();
  else if (command == GetEnabledStateCmd)
    fet->FETSimIsEnabled();
  else if (command == SetTemplateFileCmd)
    fet->SetTemplateFilename(NewValue);
  else if (command == GetTemplateFileCmd)
    fet->GetTemplateFilename();
  else if (command == SetRamoFileDirCmd)
    fet->SetRamoFileDir(NewValue);
  else if (command == GetRamoFileDirCmd)
    fet->GetRamoFileDir();
  else if (command == SetNumChanCmd)
    fet->SetNumberOfChannels(SetNumChanCmd->ConvertToInt(NewValue));
  else if (command == GetNumChanCmd)
    fet->GetNumberOfChannels();
  else if (command == SetTimeBinCmd)
    fet->SetTimeBins(SetTimeBinCmd->ConvertToDimensionedDouble(NewValue));
  else if (command == GetTimeBinCmd)
    fet->GetTimeBins();
  else if (command == SetDecayTimeCmd)
    fet->SetDecayTime(SetDecayTimeCmd->ConvertToDimensionedDouble(NewValue));
  else if (command == GetDecayTimeCmd)
    fet->GetDecayTime();
  else if (command == SetUnitTimeCmd)
    fet->SetUnitTime(SetUnitTimeCmd->ConvertToDimensionedDouble(NewValue));
  else if (command == GetUnitTimeCmd)
    fet->GetUnitTime();
  else if (command == SetPreTrigCmd)
    fet->SetPreTrig(SetPreTrigCmd->ConvertToDimensionedDouble(NewValue));
  else if (command == GetPreTrigCmd)
    fet->GetPreTrig();
  else if (command == UpdateCmd)
    fet->Build();
}

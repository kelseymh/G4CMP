/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#include "ChargeFETDigitizerModule.hh"
#include "ChargeFETDigitizerMessenger.hh"
#include "G4CMPElectrodeHit.hh"
#include "G4CMPMeshElectricField.hh"
#include "G4SystemOfUnits.hh"
#include "G4VDigitizerModule.hh"
#include "G4String.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include <sstream>

ChargeFETDigitizerModule::ChargeFETDigitizerModule(G4String modName) :
  G4VDigitizerModule(modName), messenger(new ChargeFETDigitizerMessenger(this)),
  decayTime(40e-6*s), dt(800e-9*s), preTrig(4096e-7*s), numChannels(4),
  timeBins(4096), enabledForSD(false), rereadConfigFile(true),
  rebuildFETTemplates(true), rebuildRamoFields(true),
  outputFilename("FETOutput"),
  configFilename("config/G4CMP/FETSim/ConstantsFET"),
  templateFilename("config/G4CMP/FETSim/FETTemplates"),
  ramoFileDir("config/G4CMP/FETSim")
{}

ChargeFETDigitizerModule::ChargeFETDigitizerModule() :
  G4VDigitizerModule("NoSim"), messenger(nullptr),
  decayTime(40e-6*s), dt(800e-9*s), preTrig(4096e-7*s), numChannels(4),
  timeBins(4096), enabledForSD(false), rereadConfigFile(true),
  rebuildFETTemplates(true), rebuildRamoFields(true),
  outputFilename("FETOutput"),
  configFilename("config/G4CMP/FETSim/ConstantsFET"),
  templateFilename("config/G4CMP/FETSim/FETTemplates"),
  ramoFileDir("config/G4CMP/FETSim")
{}

ChargeFETDigitizerModule::~ChargeFETDigitizerModule()
{
  delete messenger;
  if (outputFile.is_open()) outputFile.close();
  if (!outputFile.good()) {
    G4ExceptionDescription msg;
    msg << "Error closing output file, " << outputFilename << ".\n"
        << "Expect bad things like loss of data.";
    G4Exception("ChargeFETDigitizerModule::~ChargeFETDigitizerModule",
                "Charge005", FatalException, msg);
  }
}

void ChargeFETDigitizerModule::Build()
{
  SetOutputFile(outputFilename);
  if (rereadConfigFile)
    ReadFETConstantsFile();
  if (rebuildFETTemplates)
    BuildFETTemplates();
  if (rebuildRamoFields)
    BuildRamoFields();
}

void ChargeFETDigitizerModule::Digitize()
{
  if (!enabledForSD) return;
  G4HCofThisEvent* HCE =
    G4RunManager::GetRunManager()->GetCurrentEvent()->GetHCofThisEvent();
  G4SDManager* fSDM = G4SDManager::GetSDMpointer();
  G4int HCID = fSDM->GetCollectionID("G4CMPElectrodeHit");
  G4CMPElectrodeHitsCollection* hitCol =
    static_cast<G4CMPElectrodeHitsCollection*>(HCE->GetHC(HCID));
  vector<G4CMPElectrodeHit*>* hitVec = hitCol->GetVector();

  vector<G4double> scaleFactors(numChannels,0);
  G4double position[4] = {0.,0.,0.,0.};
  G4ThreeVector vecPosition;
  for(size_t chan = 0; chan < numChannels; ++chan) {
    for(size_t hitIdx=0; hitIdx < hitVec->size(); ++hitIdx) {
      if(hitVec->at(hitIdx)->GetParticleName()=="G4CMPDriftElectron" ||
         hitVec->at(hitIdx)->GetParticleName()=="G4CMPDriftHole") {
        vecPosition = hitVec->at(hitIdx)->GetFinalPosition();
        position[0] = vecPosition.getX();
        position[1] = vecPosition.getY();
        position[2] = vecPosition.getZ();
        if(hitVec->at(hitIdx)->GetParticleName()=="G4CMPDriftElectron")
          scaleFactors[chan] -= -1 * RamoFields[chan].GetPotential(position);
        else
          scaleFactors[chan] -= 1 * RamoFields[chan].GetPotential(position);
      }
    }
  }

  vector<vector<G4double> > FETTraces(CalculateTraces(scaleFactors));
  G4int runID = G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID();
  G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  WriteFETTraces(FETTraces, runID, eventID);
}

void ChargeFETDigitizerModule::PostProcess(const G4String& fileName)
{
  std::ifstream input(fileName);
  if (!input.good()) {
    G4ExceptionDescription msg;
    msg << "Error reading data input file from " << fileName;
    G4Exception("ChargeFETDigitizerModule::PostProcess", "Charge002",
    FatalException, msg);
  }
  G4double throw_away;
  G4String particleName;
  G4double position[4] = {0.,0.,0.,0.};
  G4double charge;
  G4int RunID, EventID;
  vector<G4double> scaleFactors(numChannels,0);

  G4String line;
  G4String entry;
  std::getline(input, line); //Grab column headers first
  while (!std::getline(input, line).eof()) {
    std::istringstream ssLine(line);

    std::getline(ssLine,entry,',');
    std::istringstream(entry) >> RunID;

    std::getline(ssLine,entry,',');
    std::istringstream(entry) >> EventID;

    std::getline(ssLine,entry,',');
    std::istringstream(entry) >> throw_away;

    std::getline(ssLine,entry,',');
    std::istringstream(entry) >> particleName;

    for (size_t i=0; i<6; ++i) {
      std::getline(ssLine,entry,',');
      std::istringstream(entry) >> throw_away;
    }

    for (size_t i=0; i<3; ++i) {
      std::getline(ssLine,entry,',');
      std::istringstream(entry) >> position[i];
      position[i] *= m;
    }

    std::getline(ssLine,entry,',');
    std::istringstream(entry) >> throw_away;

    if (particleName == "G4CMPDriftElectron") {
      charge = -1;
    } else if (particleName == "G4CMPDriftElectron") {
      charge = 1;
    } else {
      continue;
    }
    for(size_t chan = 0; chan < numChannels; ++chan) {
      scaleFactors[chan] -= charge*RamoFields[chan].GetPotential(position);
    }
  }

  vector<vector<G4double> > FETTraces(CalculateTraces(scaleFactors));
  WriteFETTraces(FETTraces, RunID, EventID);
}

vector<vector<G4double> > ChargeFETDigitizerModule::CalculateTraces(
                                    const vector<G4double>& scaleFactors)
{
  vector<vector<G4double> > FETTraces(numChannels,vector<G4double>(timeBins,0.));
  for(size_t chan=0; chan < numChannels; ++chan)
    for(size_t cross=0; cross < numChannels; ++cross)
      for(size_t bin=0; bin < timeBins; ++bin)
        FETTraces[chan][bin] += scaleFactors[cross]*FETTemplates[chan][cross][bin];
  return FETTraces;
}

void ChargeFETDigitizerModule::ReadFETConstantsFile()
{
  constantsFile.open(configFilename);
  if (!constantsFile.good()) {
    G4ExceptionDescription msg;
    msg << "Error reading FET constants file from " << configFilename;
    G4Exception("ChargeFETDigitizerModule::ReadFETConstantsFile", "Charge001",
		FatalException, msg);
  }
  G4String buffer;
  G4String varName;
  G4String varVal;
  G4bool comment = false;
  size_t pos;

  while(!constantsFile.eof()) {
    getline(constantsFile, buffer);
    if((!comment) && (buffer.find("//") != std::string::npos)) {
      if(buffer.find("//") == 0)
        continue;
      else
        buffer = buffer.substr(0, buffer.find("//"));
    } else if((!comment) && (buffer.find("/*") != std::string::npos)) {
      buffer.erase(buffer.find("/*"), buffer.find("*/") + 2);
      if(buffer.find("*/") == std::string::npos)
        comment = true;
    } else if(buffer.find("*/") != std::string::npos && comment) {
      comment = false;
      buffer.erase(buffer.find("*/"), 2);
    }

    if(!comment && !buffer.empty()) {
      while(buffer.find(' ') != std::string::npos) //Erase all spaces
        buffer.erase(buffer.find(' '), 1);

      pos = buffer.find('='); //To split buffer into lhs and rhs
      if(pos == std::string::npos)
        continue;
      else {
        varName = buffer.substr(0, pos++);
        varVal = buffer.substr(pos);

        if(*(--varVal.end()) == ';') //Strip ; from end
          varVal.erase(--varVal.end());

        if(varName == "numChannels")
          SetNumberOfChannels(atoi(varVal));
        else if(varName == "timeBins")
          SetTimeBins(atoi(varVal));
        else if(varName == "decayTime")
          SetDecayTime(atof(varVal)*s);
        else if(varName == "dt")
          SetUnitTime(atof(varVal)*s);
        else if(varName == "preTrig")
          SetPreTrig(atof(varVal)*s);
        else if(varName == "templateFilename")
          SetTemplateFilename(varVal.substr(1,varVal.length()-2)); //strip quotes
        else if(varName == "ramoFileDir") {
          varVal = varVal.substr(1,varVal.length()-2); //strip quotes
          if (*(--varVal.end()) == '/') varVal.erase(--varVal.end()); //strip trailing slash
          SetRamoFileDir(varVal);
        }
        else
          G4cout << "FETSim Warning: Variable " << varName
                 << " is unused." << G4endl;
      }
    }
  }
  constantsFile.close();
  rereadConfigFile = false;
}

void ChargeFETDigitizerModule::BuildFETTemplates()
{
  FETTemplates = vector<vector<vector<G4double> > >
    (numChannels, vector<vector<G4double> >
    (numChannels, vector<G4double>(timeBins, 0) ) );
  templateFile.open(templateFilename.c_str());
  if(templateFile.good()) {
    for(size_t i=0; i<numChannels; ++i)
      for(size_t j=0; j<numChannels; ++j)
        for(size_t k=0; k<timeBins; ++k)
          templateFile >> FETTemplates[i][j][k];
  } else {
    G4Exception("ChargeFETDigitizerModule::BuildFETTemplate", "Charge007",
		JustWarning,
	"Reading from template file failed. Using default pulse templates.");

    for(size_t i=0; i<numChannels; ++i) {
      size_t ndt = static_cast<size_t>(preTrig/dt);
      for(size_t j=0; j<ndt; ++j)
        FETTemplates[i][i][j] = 0;
      for(size_t k=1; k<timeBins-ndt+1; ++k)
        FETTemplates[i][i][k+ndt-1] = exp(-k*dt/decayTime);
    }
  }
  templateFile.close();
  rebuildFETTemplates = false;
}

void ChargeFETDigitizerModule::BuildRamoFields()
{
  if (RamoFields.size()) RamoFields.clear();

  for(size_t i=0; i < numChannels; ++i) {
    std::stringstream name;
    name << ramoFileDir << "/EpotRamoChan" << i+1;
    std::ifstream ramoFile(name.str().c_str());
    if(ramoFile.good()) {
      ramoFile.close();
      RamoFields.emplace_back(name.str());
    } else {
      ramoFile.close();
      G4cerr << "ChargeFETDigitizerModule::BuildRamoFields(): ERROR: Could"
        << " not open Ramo files for each FET channel." << G4endl;
    }
  }
  rebuildRamoFields = false;
}

void ChargeFETDigitizerModule::WriteFETTraces(
  const vector<vector<G4double> >& traces, G4int RunID, G4int EventID)
{
  for(size_t chan = 0; chan < numChannels; ++chan) {
    outputFile << RunID << "," << EventID << "," << chan+1 << ",";
    for(size_t bin = 0; bin < timeBins-1; ++bin) {
      outputFile << traces[chan][bin] << ",";
    }
    outputFile << traces[chan][timeBins-1] << "\n";
  }
}

void ChargeFETDigitizerModule::EnableFETSim()
{
  enabledForSD = true;
  if (RamoFields.size() == 0) { // Need to initiate first build.
    Build();
  }
}

void ChargeFETDigitizerModule::SetOutputFile(const G4String& fn)
{
  if (outputFilename != fn) {
    if (outputFile.is_open()) outputFile.close();
    outputFilename = fn;
    outputFile.open(outputFilename, std::ios_base::app);
    if (!outputFile.good()) {
      G4ExceptionDescription msg;
      msg << "Error opening output file, " << outputFilename << ".\n"
          << "Will continue simulation.";
      G4Exception("ChargeFETDigitizerModule::SetOutputFile", "Charge006",
                  JustWarning, msg);
      outputFile.close();
    } else {
      outputFile << "Run ID,Event ID,Channel,Pulse (4096 bins)" << G4endl;
    }
  }
}

void ChargeFETDigitizerModule::SetConfigFilename(const G4String& name)
{
  if (configFilename == name) return;
  configFilename = name;
  rereadConfigFile = true;
}

void ChargeFETDigitizerModule::SetTemplateFilename(const G4String& name)
{
  if (templateFilename == name) return;
  templateFilename = name;
  rebuildFETTemplates = true;
}

void ChargeFETDigitizerModule::SetRamoFileDir(const G4String& name)
{
  if (ramoFileDir == name) return;
  ramoFileDir = name;
  rebuildRamoFields = true;
}

void ChargeFETDigitizerModule::SetNumberOfChannels(size_t n)
{
  if (numChannels == n) return;
  numChannels = n;
  rebuildRamoFields = true;
  rebuildFETTemplates = true;
}

void ChargeFETDigitizerModule::SetTimeBins(size_t n)
{
  if (timeBins == n) return;
  timeBins = n;
  rebuildFETTemplates = true;
}

void ChargeFETDigitizerModule::SetDecayTime(G4double n)
{
  if (decayTime == n) return;
  decayTime = n;
  rebuildFETTemplates = true;
}

void ChargeFETDigitizerModule::SetUnitTime(G4double n)
{
  if (dt == n) return;
  dt = n;
  rebuildFETTemplates = true;
}

void ChargeFETDigitizerModule::SetPreTrig(G4double n)
{
  if (preTrig == n) return;
  preTrig = n;
  rebuildFETTemplates = true;
}

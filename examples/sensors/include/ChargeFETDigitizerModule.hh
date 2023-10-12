/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef CHARGEFETDIGITIZERMODULE_HH
#define CHARGEFETDIGITIZERMODULE_HH

#include "G4VDigitizerModule.hh"
#include <fstream>

class ChargeFETDigitizerMessenger;
class G4CMPMeshElectricField;
class G4CMPElectrodeHit;
class G4String;

using std::vector;

class ChargeFETDigitizerModule : public G4VDigitizerModule
{
  public:
    ChargeFETDigitizerModule(G4String modName);
    // Default constructor only to be used for stand-alone post-processing!
    ChargeFETDigitizerModule();
    virtual ~ChargeFETDigitizerModule();

    void Build();
    virtual void Digitize();
    void PostProcess(const G4String& fileName);

    // Methods for Messenger
    void     EnableFETSim();
    void     DisableFETSim() {enabledForSD = false;}
    G4bool   FETSimIsEnabled() const {return enabledForSD;}

    void     SetOutputFile(const G4String& name);
    G4String GetOutputFile() const {return outputFilename;}

    void     SetConfigFilename(const G4String& name);
    G4String GetConfigFilename() const {return configFilename;}

    void     SetTemplateFilename(const G4String& name);
    G4String GetTemplateFilename() const {return templateFilename;}

    void     SetRamoFileDir(const G4String& name);
    G4String GetRamoFileDir() const {return ramoFileDir;}

    void     SetNumberOfChannels(size_t n);
    size_t   GetNumberOfChannels() const {return numChannels;}

    void     SetTimeBins(size_t n);
    size_t   GetTimeBins() const {return timeBins;}

    void     SetDecayTime(G4double n);
    G4double GetDecayTime() const {return decayTime;}

    void     SetUnitTime(G4double n);
    G4double GetUnitTime() const {return dt;}

    void     SetPreTrig(G4double n);
    G4double GetPreTrig() const {return preTrig;}

  private:
    void ReadFETConstantsFile();
    void BuildFETTemplates();
    vector<vector<G4double> > CalculateTraces(const vector<G4double>& scaleFactors);
    void BuildRamoFields();
    void WriteFETTraces(const vector<vector<G4double> >& FETTraces,
                        G4int RunID, G4int EventID);

    ChargeFETDigitizerMessenger* messenger;
    // FET constants
    G4double decayTime;
    G4double dt;
    G4double preTrig;
    size_t numChannels;
    size_t timeBins;
    // Enable/Disable FETSim during sim
    G4bool enabledForSD;
    // Internal flags to not waste time on unnecessary recalculating
    G4bool rereadConfigFile;
    G4bool rebuildFETTemplates;
    G4bool rebuildRamoFields;
    // File Stuff
    std::ofstream outputFile;
    std::ifstream constantsFile;
    std::ifstream templateFile;
    G4String outputFilename;
    G4String configFilename;
    G4String templateFilename;
    G4String ramoFileDir;
    // FETSim Quantities
    vector<vector<vector<G4double> > > FETTemplates; //4x4x4096 = 4 channels w/ cross-talk terms
    vector<G4CMPMeshElectricField> RamoFields;
};

#endif // CHARGEFETDIGITIZERMODULE_HH

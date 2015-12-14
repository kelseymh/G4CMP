#ifndef CHARGEFETDIGITIZERMODULE_HH
#define CHARGEFETDIGITIZERMODULE_HH

#include "G4VDigitizerModule.hh"
#include <vector>
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
    // Default constructor only to be used for post-processing!
    ChargeFETDigitizerModule();
    virtual ~ChargeFETDigitizerModule();

    void Initialize();
    virtual void Digitize();
    void PostProcess(const G4String& fileName);

    // Methods for Messenger
    void     EnableFETSim() {runLive = true;}
    void     DisableFETSim() {runLive = false;}
    G4bool   FETSimIsEnabled() const {return runLive;}

    void     SetOutputFilename(G4String name) {outputFilename = name;}
    G4String GetOutputFilename() const {return outputFilename;}

    void     SetConfigFilename(G4String name) {configFilename = name;}
    G4String GetConfigFilename() const {return configFilename;}

    void     SetTemplateFilename(G4String name) {templateFilename = name;}
    G4String GetTemplateFilename() const {return templateFilename;}

    void     SetRamoFileDir(G4String name) {ramoFileDir = name;}
    G4String GetRamoFileDir() const {return ramoFileDir;}

    void     SetNumberOfChannels(G4int n) {numChannels = n;}
    G4int    GetNumberOfChannels() const {return numChannels;}

    void     SetTimeBins(G4int n) {timeBins = n;}
    G4int    GetTimeBins() const {return timeBins;}

    void     SetDecayTime(G4double n) {decayTime = n;}
    G4double GetDecayTime() const {return decayTime;}

    void     SetUnitTime(G4double n) {dt = n;}
    G4double GetUnitTime() const {return dt;}

    void     SetPreTrig(G4double n) {preTrig = n;}
    G4double GetPreTrig() const {return preTrig;}

    void     SetTemplateEnergy(G4double n) {templateEnergy = n;}
    G4double GetTemplateEnergy() const {return templateEnergy;}


  private:
    void ReadFETConstantsFile();
    void BuildFETTemplates();
    vector<vector<G4double> > CalculateTraces(const vector<G4double>& scaleFactors);
    void BuildRamoFields();
    void WriteFETTraces(const vector<vector<G4double> >& FETTraces,
                        G4int RunID, G4int EventID);

    // File Stuff
    std::ofstream outputFile;
    std::ifstream constantsFile;
    std::ifstream templateFile;
    G4String outputFilename;
    G4String configFilename;
    G4String templateFilename;
    G4String ramoFileDir;
    // FET constants
    G4double decayTime;
    G4double dt;
    G4double preTrig;
    G4double templateEnergy;
    G4int numChannels;
    G4int timeBins;
    // Enable/Disable FETSim during sim
    G4bool runLive;
    // FETSim Quantities
    vector<vector<vector<G4double> > > FETTemplates; //4x4x4096 = 4 channels w/ cross-talk terms
    vector<G4CMPMeshElectricField> RamoFields;
    ChargeFETDigitizerMessenger* messenger;
};

#endif // CHARGEFETDIGITIZERMODULE_HH

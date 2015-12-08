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
    void PostProcess(G4String fileName);

    // Methods for Messenger
    void     SetOutputFilename(G4String name) {outputFilename = name;}
    G4String GetOutputFilename() {return outputFilename;}

    void     SetTemplateFilename(G4String name) {templateFilename = name;}
    G4String GetTemplateFilename() {return templateFilename;}

    void     SetRamoFileDir(G4String name) {ramoFileDir = name;}
    G4String GetRamoFileDir() {return ramoFileDir;}

    void     SetNumberOfChannels(G4int n) {numChannels = n;}
    G4int    GetNumberOfChannels() {return numChannels;}

    void     SetTimeBins(G4int n) {timeBins = n;}
    G4int    GetTimeBins() {return timeBins;}

    void     SetDecayTime(G4double n) {decayTime = n;}
    G4double GetDecayTime() {return decayTime;}

    void     SetUnitTime(G4double n) {dt = n;}
    G4double GetUnitTime() {return dt;}

    void     SetPreTrig(G4double n) {preTrig = n;}
    G4double GetPreTrig() {return preTrig;}

    void     SetTemplateEnergy(G4double n) {templateEnergy = n;}
    G4double GetTemplateEnergy() {return templateEnergy;}


  private:
    void ReadFETConstantsFile(G4String filename);
    void BuildFETTemplates();
    vector<vector<G4double> > CalculateTraces(const vector<G4double>& scaleFactors);
    void BuildRamoFields();
    void WriteFETTraces(const vector<vector<G4double> >& FETTraces,
                        G4int RunID, G4int EventID);

    std::fstream output;
    ChargeFETDigitizerMessenger* messenger;
    G4String configFilename;
    G4String outputFilename;
    //FET constants
    G4int numChannels, timeBins;
    G4double decayTime, dt, preTrig, templateEnergy;
    G4String templateFilename;
    G4String ramoFileDir;
    //FETSim Quantities
    vector<vector<vector<G4double> > > FETTemplates; //4x4x4096 = 4 channels w/ cross-talk terms
    vector<G4CMPMeshElectricField*> RamoFields;
};

#endif // CHARGEFETDIGITIZERMODULE_HH

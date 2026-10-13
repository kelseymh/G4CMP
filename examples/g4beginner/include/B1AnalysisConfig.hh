#ifndef B1AnalysisConfig_h
#define B1AnalysisConfig_h 1

#include "globals.hh"
#include "G4String.hh"

class G4GenericMessenger;

class B1AnalysisConfig
{
  public:
    static B1AnalysisConfig* Instance();
    ~B1AnalysisConfig();

    void SetSpecies(const G4String& species);
    const G4String& GetSpecies() const;

    G4GenericMessenger* GetMessenger() const { return Instance()->fMessenger; }

    G4bool IsElectronMode() const;
    G4bool IsPhononMode() const;

    G4String GetOutputFileName() const;
  private:
    B1AnalysisConfig();
    B1AnalysisConfig(const B1AnalysisConfig&);  // To clone from master

    G4String fSpecies;
    G4GenericMessenger* fMessenger;
};

#endif

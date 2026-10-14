#include "B1AnalysisConfig.hh"

#include "G4GenericMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"
#include "G4ios.hh"

B1AnalysisConfig* B1AnalysisConfig::Instance()
{
  // ripped straight from G4CMPConfigManager
  static const B1AnalysisConfig* masterInstance = 0;
  static G4ThreadLocal B1AnalysisConfig* theInstance = 0;

  if (!theInstance) {
    if (!G4Threading::IsWorkerThread()) { // Master or sequential
      theInstance = new B1AnalysisConfig;
      masterInstance = theInstance;
    } else {          // Workers copy from master
      theInstance = new B1AnalysisConfig(*masterInstance);
    }
  }

  return theInstance;
}

B1AnalysisConfig::B1AnalysisConfig()
: fSpecies("electrons"),
  fMessenger(nullptr)
{
  fMessenger = new G4GenericMessenger(this, "/analysis/", "Analysis control");

  auto& speciesCmd =
    fMessenger->DeclareMethod("setSpecies", &B1AnalysisConfig::SetSpecies,
                              "Set species to analyze: electrons or phonons");
  speciesCmd.SetParameterName("species", false);
}

B1AnalysisConfig::B1AnalysisConfig(const B1AnalysisConfig& rhs)
: fSpecies(rhs.GetSpecies()),
  fMessenger(rhs.GetMessenger())
{
  if (fMessenger) {
  auto& speciesCmd =
      fMessenger->DeclareMethod("setSpecies", &B1AnalysisConfig::SetSpecies,
                                "Set species to analyze: electrons or phonons");
    speciesCmd.SetParameterName("species", false);
  }
}

B1AnalysisConfig::~B1AnalysisConfig()
{
  delete fMessenger;
  fMessenger = 0;
}

void B1AnalysisConfig::SetSpecies(const G4String& species)
{
  if (species == "electrons" || species == "phonons") {
    Instance()->fSpecies = species;
    G4cout << "Analysis species set to: " << fSpecies << G4endl;
  } else {
    G4cout << "Unknown species mode: " << species
           << ". Valid options are: electrons, phonons" << G4endl;
  }
}

const G4String& B1AnalysisConfig::GetSpecies() const
{
  return Instance()->fSpecies;
}

G4bool B1AnalysisConfig::IsElectronMode() const
{
  return Instance()->fSpecies == "electrons";
}

G4bool B1AnalysisConfig::IsPhononMode() const
{
  return Instance()->fSpecies == "phonons";
}

G4String B1AnalysisConfig::GetOutputFileName() const
{
  if (Instance()->fSpecies == "electrons") return "g4cmp_electrons.root";
  if (Instance()->fSpecies == "phonons")   return "g4cmp_phonons.root";
  return "g4cmp_output.root";
}

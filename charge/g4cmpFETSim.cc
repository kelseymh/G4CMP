#include "ChargeFETDigitizerModule.hh"
#include "G4SystemOfUnits.hh"

int main(int argc, char** argv) {
  G4String filename;
  if (argc == 1) {
    G4cout << "Enter path to data file to be processed: " << G4endl;
    G4cin >> filename;
  } else {
    filename = argv[1];
  }

  ChargeFETDigitizerModule fetsim;
  if (argc > 2)
    fetsim.SetOutputFilename(argv[2]);
  fetsim.Initialize();
  fetsim.PostProcess(filename);

  return 0;
}

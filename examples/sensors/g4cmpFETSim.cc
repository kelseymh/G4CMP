/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#include "ChargeFETDigitizerModule.hh"

int main(int argc, char** argv) {
  G4String filename;
  if (argc == 1) {
    G4cout << "Enter path to data file to be processed: " << G4endl;
    G4cin >> filename;
  } else {
    filename = argv[1];
  }

  ChargeFETDigitizerModule fetsim;
  if (argc > 2) {
    fetsim.SetOutputFile(argv[2]);
  } else {
    fetsim.SetOutputFile("FETOutput");
  }

  fetsim.Build();
  fetsim.PostProcess(filename);

  return 0;
}

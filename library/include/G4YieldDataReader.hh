/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4YieldDataReader.hh
/// \Reads the yield data for Sarkis model from /CrystalMaps/NIEL/SarkisSiIonYield.txt, populates the relevant data and store it in a other file to be used in G4CMPSarkisNIEL. This could be extended for other models that should be implemented the same way. 
///
// 20230814  David Sadek


#ifndef G4YIELDDATAREADER_HH
#define G4YIELDDATAREADER_HH

#include "G4String.hh"
#include "G4PhysicsLogVector.hh"
#include "globals.hh"
#include <vector>

class G4YieldDataReader {
public:
    G4YieldDataReader(const G4String& dataFileName, const G4String& outputFileName);
    ~G4YieldDataReader();

    void ReadDataAndStore();

private:
    G4String fDataFileName;
    G4String fOutputFileName;
};

#endif // G4YIELDDATAREADER_HH


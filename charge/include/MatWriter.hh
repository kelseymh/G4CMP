/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef MatWriter_h
#define MatWriter_h

#include "globals.hh"
#include <vector>
using std::vector;

class MatWriter
{
    public:
       MatWriter(G4String outFile, G4String varName, G4int rows, G4int columns, vector<G4double> data, G4bool append);
       ~MatWriter() {};
};
#endif

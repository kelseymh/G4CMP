/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id: 4de7ad9cca0c224be0c33719a95cbeeb647113af $
//
#ifndef XLogicalAtomicLatticeDiamond_h
#define XLogicalAtomicLatticeDiamond_h

#include <iostream>
#include <fstream>
#include <string>
#include "G4ThreeVector.hh"
#include "XLogicalAtomicLattice.hh"

using namespace std;

class XLogicalAtomicLatticeDiamond: public XLogicalAtomicLattice
{

private:
    void InitializeXLogicalAtomicLatticeDiamond();

public:    
    XLogicalAtomicLatticeDiamond();
    ~XLogicalAtomicLatticeDiamond();
};

#endif

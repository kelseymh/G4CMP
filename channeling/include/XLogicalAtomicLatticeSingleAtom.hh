/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id: f01783235cfdccb92049cc97296a25e1d1fe1777 $
//
#ifndef XLogicalAtomicLatticeSingleAtom_h
#define XLogicalAtomicLatticeSingleAtom_h

#include <iostream>
#include <fstream>
#include <string>
#include "G4ThreeVector.hh"
#include "XLogicalAtomicLattice.hh"

using namespace std;

class XLogicalAtomicLatticeSingleAtom: public XLogicalAtomicLattice
{

private:
    void InitializeXLogicalAtomicLatticeSingleAtom();

public:    
    XLogicalAtomicLatticeSingleAtom();
    ~XLogicalAtomicLatticeSingleAtom();
};

#endif

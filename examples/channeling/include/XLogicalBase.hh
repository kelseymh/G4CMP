/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id: acc1dca33b925e511839f9f46b51efb19e784f67 $
//
#ifndef XLogicalBase_h
#define XLogicalBase_h

#include <iostream>
#include <fstream>
#include <string>
#include "G4ThreeVector.hh"
#include "XLogicalAtomicLattice.hh"
#include "G4Element.hh"

using namespace std;

class XLogicalBase{

private:
    XLogicalAtomicLattice* fLattice;
    G4Element* fElement;
    
public:
    // Retrieval methods
    XLogicalAtomicLattice* GetLattice();
    G4Element* GetElement();
    
    // Set methods
    void SetLattice(XLogicalAtomicLattice*);
    void SetElement(G4Element*);
    
    // Calculation methods
    // ints == Miller indexes
    virtual G4double ComputeAtomicFormFactor(); //Kittel - chapter 2 Eq. (42) for single atomic kind
    G4complex ComputeStructureFactorSingleAtomicKind(G4int,G4int,G4int);  //Kittel - chapter 2 Eq. (46) for single atomic kind
  
    XLogicalBase(G4Element*,XLogicalAtomicLattice*);
    XLogicalBase();
    ~XLogicalBase();
};

#endif

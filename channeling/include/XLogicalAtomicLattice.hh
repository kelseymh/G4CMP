/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id: 1c81d6ffd84ba064209a4d3b4c84806e0ca7889d $
//
#ifndef XLogicalAtomicLattice_h
#define XLogicalAtomicLattice_h

#include <iostream>
#include <fstream>
#include <string>
#include "G4ThreeVector.hh"

#ifndef MAXLATTICEATOMS
#define MAXLATTICEATOMS 64
#endif

using namespace std;

class XLogicalAtomicLattice{

private:
    // position of the atoms are saved in unit cell system, i.e MIN 0. & MAX 1.
    G4ThreeVector fLatticeAtomPosition[MAXLATTICEATOMS];
    G4int fLatticeAtomNumber;
    
public:    
    void InitializeXLogicalAtomicLattice();

    // Get methods
    G4ThreeVector GetAtomPosition(G4int i);
    G4int GetLatticeNumberOfAtoms();
    
    // Set methods
    void AddAtom(G4ThreeVector);
    void DeleteAtom(G4ThreeVector);
    

    // Calculation methods
    // ints == Miller indexes
    G4complex ComputeGeometricalStructureFactorSingleKind(G4int,G4int,G4int);

    // Definition methods
    XLogicalAtomicLattice();
    ~XLogicalAtomicLattice();
};

#endif

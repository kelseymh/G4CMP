/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#include "XLogicalAtomicLatticeDiamond.hh"
#include "G4PhysicalConstants.hh"
#include <cmath>

XLogicalAtomicLatticeDiamond::XLogicalAtomicLatticeDiamond(){
    InitializeXLogicalAtomicLatticeDiamond();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XLogicalAtomicLatticeDiamond::~XLogicalAtomicLatticeDiamond(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XLogicalAtomicLatticeDiamond::InitializeXLogicalAtomicLatticeDiamond(){
    InitializeXLogicalAtomicLattice();

    for(unsigned int i=0;i<2;i++)
    {
        AddAtom(G4ThreeVector(0.0+0.25*i,0.0+0.25*i,0.0+0.25*i));
        AddAtom(G4ThreeVector(0.5+0.25*i,0.5+0.25*i,0.0+0.25*i));
        AddAtom(G4ThreeVector(0.0+0.25*i,0.5+0.25*i,0.5+0.25*i));
        AddAtom(G4ThreeVector(0.5+0.25*i,0.0+0.25*i,0.5+0.25*i));
    }
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

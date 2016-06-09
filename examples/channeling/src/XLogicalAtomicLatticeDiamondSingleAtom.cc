/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#include "XLogicalAtomicLatticeSingleAtom.hh"
#include "G4PhysicalConstants.hh"
#include <cmath>

XLogicalAtomicLatticeSingleAtom::XLogicalAtomicLatticeSingleAtom(){
    InitializeXLogicalAtomicLatticeSingleAtom();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XLogicalAtomicLatticeSingleAtom::~XLogicalAtomicLatticeSingleAtom(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XLogicalAtomicLatticeSingleAtom::InitializeXLogicalAtomicLatticeSingleAtom(){
    InitializeXLogicalAtomicLattice();
    AddAtom(G4ThreeVector(0.,0.,0.));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

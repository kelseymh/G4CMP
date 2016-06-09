/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#include "XLogicalAtomicLattice.hh"
#include "G4PhysicalConstants.hh"
#include <cmath>

XLogicalAtomicLattice::XLogicalAtomicLattice(){
    InitializeXLogicalAtomicLattice();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XLogicalAtomicLattice::~XLogicalAtomicLattice(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XLogicalAtomicLattice::InitializeXLogicalAtomicLattice(){
    fLatticeAtomNumber = 1;
    for(G4int i=0;i<MAXLATTICEATOMS;i++) fLatticeAtomPosition[i] = G4ThreeVector(0.,0.,0.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector XLogicalAtomicLattice::GetAtomPosition(G4int i){
    if(i<fLatticeAtomNumber){
        return fLatticeAtomPosition[i];
    }
    else{
        G4cout << "XLogicalAtomicLattice::GetAtomPosition - atom " << i << " does not exist!!" <<std::endl;
    }
    return G4ThreeVector(-1.,-1.,-1.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int XLogicalAtomicLattice::GetLatticeNumberOfAtoms(){
    return fLatticeAtomNumber;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XLogicalAtomicLattice::AddAtom(G4ThreeVector vAtomPosition){
    fLatticeAtomNumber++;
    //Add an atom to the lattice
    fLatticeAtomPosition[fLatticeAtomNumber - 1] = vAtomPosition;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XLogicalAtomicLattice::DeleteAtom(G4ThreeVector vAtomPosition){
    //Delete atoms in the lattice in the selected position
    
    for(G4int i=0;i<fLatticeAtomNumber;i++)
        if(vAtomPosition == fLatticeAtomPosition[i])
        {
            for(G4int j=(i+1);j<fLatticeAtomNumber;j++)
            {
                fLatticeAtomPosition[j-1]=fLatticeAtomPosition[j];
            }
            i--;
            fLatticeAtomNumber--;
        }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4complex XLogicalAtomicLattice::ComputeGeometricalStructureFactorSingleKind(G4int h,G4int k ,G4int l){
    G4double vTempDouble = 0.;
    G4complex vResult = G4complex(0.,0.);

    for(G4int i=0;i<fLatticeAtomNumber;i++)
    {
        vTempDouble = 0.0;
        vTempDouble += h * fLatticeAtomPosition[i].x();
        vTempDouble += k * fLatticeAtomPosition[i].y();
        vTempDouble += l * fLatticeAtomPosition[i].z();
        vResult += G4complex(cos(2 * M_PI * vTempDouble),sin(2 * M_PI * vTempDouble));
    }

    return vResult;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

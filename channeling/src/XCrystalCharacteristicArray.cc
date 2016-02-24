/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#include "XCrystalCharacteristicArray.hh"

XCrystalCharacteristicArray::XCrystalCharacteristicArray(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XCrystalCharacteristicArray::~XCrystalCharacteristicArray(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector XCrystalCharacteristicArray::ComputeEC(G4ThreeVector vPosition,XPhysicalLattice* vLattice){
    
    G4ThreeVector vValue = G4ThreeVector(0.,0.,0.);
    
    if(fCharacteristicVector.size()!=0){
        for(unsigned int i=0;i<fCharacteristicVector.size();i++){
            vValue += fCharacteristicVector.at(i)->GetEC(vPosition,vLattice);
        }
        return (vValue * G4double(1./G4double(fCharacteristicVector.size())));
    }
    
    return G4ThreeVector(DBL_MAX,DBL_MAX,DBL_MAX);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector XCrystalCharacteristicArray::ComputePositionInUnitCell(G4ThreeVector vPosition,XPhysicalLattice*vLattice){
    
    G4ThreeVector vPositionInTheCellPrevious = G4ThreeVector(0.,0.,0.);
    G4ThreeVector vPositionInTheCell = G4ThreeVector(0.,0.,0.);
    
    if(fCharacteristicVector.size()!=0){
        for(unsigned int i=0;i<fCharacteristicVector.size();i++){
            if(vPositionInTheCellPrevious != vPositionInTheCell){
                return G4ThreeVector(DBL_MAX,DBL_MAX,DBL_MAX);
            }
            vPositionInTheCellPrevious = vPositionInTheCell;
            vPositionInTheCell = fCharacteristicVector.at(i)->ComputePositionInUnitCell(vPosition,vLattice);
        }
        return vPositionInTheCell;
    }
    
    return G4ThreeVector(DBL_MAX,DBL_MAX,DBL_MAX);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


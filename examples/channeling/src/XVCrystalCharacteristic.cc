/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#include "XVCrystalCharacteristic.hh"

XVCrystalCharacteristic::XVCrystalCharacteristic(){
    fLatticeManager = XLatticeManager3::GetXLatticeManager();
    
    fMaximum = DBL_MAX;
    fMinimum = DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XVCrystalCharacteristic::~XVCrystalCharacteristic(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XPhysicalLattice* XVCrystalCharacteristic::GetXPhysicalLattice(G4VPhysicalVolume* vVolume){
    return fLatticeManager->GetXPhysicalLattice(vVolume);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XUnitCell* XVCrystalCharacteristic::GetXUnitCell(G4VPhysicalVolume* vVolume){
    return GetXPhysicalLattice(vVolume)->GetXUnitCell();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XLogicalLattice* XVCrystalCharacteristic::GetLogicalLattice(G4VPhysicalVolume* vVolume){
    return GetXPhysicalLattice(vVolume)->GetLogicalLattice();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XVCrystalCharacteristic::InitializePhysicalLattice(XPhysicalLattice* vLattice){
    if(fPhysicalLattice != vLattice){
        fPhysicalLattice = vLattice;
        InitializeVector();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector XVCrystalCharacteristic::GetEC(G4ThreeVector vPosition,XPhysicalLattice* vLattice){
    if(IsInitialized(vLattice)){
        return ComputeECFromVector(vPosition);
    }
    else{
        return ComputeEC(vPosition,vLattice);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XVCrystalCharacteristic::ComputeTFScreeningRadius(XPhysicalLattice* vLattice){
    
    G4double vTFSR = CLHEP::Bohr_radius * 0.88534;
    
    vTFSR /= (std::pow(vLattice->GetXUnitCell()->GetBase(0)->GetElement()->GetZ(),0.333333333));
    
    //    if(vLattice->GetParticleDefinition()->GetParticleName() == "proton"){
    //        vTFSR /= (std::pow(vLattice->GetMaterial()->GetZ(),0.333333333));
    //    }
    //    else{
    //        vTFSR /= (std::pow(vLattice->GetMaterial()->GetZ(),0.23) + std::pow(vLattice->GetParticleDefinition()->GetPDGCharge(),0.23));
    //    }
    
    return vTFSR;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector XVCrystalCharacteristic::ComputePositionInUnitCell(G4ThreeVector,XPhysicalLattice*){
    return G4ThreeVector(-1.,-1.,-1.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XVCrystalCharacteristic::GetMaximum(XPhysicalLattice* vLattice){
    if(fMaximum == DBL_MAX){
        fMaximum = ComputeMaximum(vLattice);
    }
    return fMaximum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XVCrystalCharacteristic::GetMinimum(XPhysicalLattice* vLattice){
    if(fMinimum == DBL_MAX){
        fMinimum = ComputeMinimum(vLattice);
    }
    return fMinimum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XVCrystalCharacteristic::ComputeMaximum(XPhysicalLattice* vLattice){
    unsigned int vPrecisionX = 1024;
    unsigned int vPrecisionY = 1024;
    unsigned int vPrecisionZ = 1024;
    
    G4VPhysicalVolume* vVolume = vLattice->GetVolume();
    G4double vStepX = GetXUnitCell(vVolume)->GetSize().x() / vPrecisionX;
    G4double vStepY = GetXUnitCell(vVolume)->GetSize().y() / vPrecisionY;
    G4double vStepZ = GetXUnitCell(vVolume)->GetSize().z() / vPrecisionZ;
    
    G4double vMaximum = -DBL_MAX;
    G4double vValue;
    
    for(unsigned int i=0;i<vPrecisionX;i++){
        for(unsigned int j=0;j<vPrecisionY;j++){
            for(unsigned int k=0;k<vPrecisionZ;k++){
                if( (vValue = GetEC(G4ThreeVector(vStepX * i,vStepY * i,vStepZ * i),vLattice).mag() ) > vMaximum) vMaximum = vValue;
            }
        }
    }
    
    return vMaximum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XVCrystalCharacteristic::ComputeMinimum(XPhysicalLattice* vLattice){
    unsigned int vPrecisionX = 1024;
    unsigned int vPrecisionY = 1024;
    unsigned int vPrecisionZ = 1024;
    
    G4VPhysicalVolume* vVolume = vLattice->GetVolume();
    G4double vStepX = GetXUnitCell(vVolume)->GetSize().x() / vPrecisionX;
    G4double vStepY = GetXUnitCell(vVolume)->GetSize().y() / vPrecisionY;
    G4double vStepZ = GetXUnitCell(vVolume)->GetSize().z() / vPrecisionZ;
    
    G4double vMinimum = DBL_MAX;
    G4double vValue;
    
    for(unsigned int i=0;i<vPrecisionX;i++){
        for(unsigned int j=0;j<vPrecisionY;j++){
            for(unsigned int k=0;k<vPrecisionZ;k++){
                if( (vValue = GetEC(G4ThreeVector(vStepX * i,vStepY * i,vStepZ * i),vLattice).mag() ) < vMinimum) vMinimum = vValue;
            }
        }
    }
    
    return vMinimum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool XVCrystalCharacteristic::IsInitialized(XPhysicalLattice* vLattice){
    if(vLattice == fPhysicalLattice){
        return true;
    }
    else{
        return false;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

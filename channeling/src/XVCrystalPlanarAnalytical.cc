/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#include "XVCrystalPlanarAnalytical.hh"
#include "G4PhysicsLinearVector.hh"

XVCrystalPlanarAnalytical::XVCrystalPlanarAnalytical(){
    fNumberOfPlanes = 4;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XVCrystalPlanarAnalytical::~XVCrystalPlanarAnalytical(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XVCrystalPlanarAnalytical::SetNumberOfPlanes(G4int vNumberOfPlanes){
    fNumberOfPlanes = vNumberOfPlanes;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int XVCrystalPlanarAnalytical::GetNumberOfPlanes(){
    return fNumberOfPlanes;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector XVCrystalPlanarAnalytical::ComputeEC(G4ThreeVector vPositionVector,XPhysicalLattice* vLattice){
    G4VPhysicalVolume* vVolume = vLattice->GetVolume();
    
    G4double vInterplanarDistance = GetXPhysicalLattice(vVolume)->ComputeInterplanarPeriod();
    
    G4double vPosition = ComputePositionInUnitCell(vPositionVector,vLattice).x();
    
    G4double vValue = 0.;
    for(int i=-int(GetNumberOfPlanes()/2);i<=+int(GetNumberOfPlanes()/2);i++){
        vValue += ComputeECForSinglePlane( ( vPosition + G4double(i) ) * vInterplanarDistance , vLattice );
    }
    
    return G4ThreeVector(vValue,0.,0.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector XVCrystalPlanarAnalytical::ComputeECFromVector(G4ThreeVector vPosition){
    G4double vInterplanarPeriod = fPhysicalLattice->ComputeInterplanarPeriod();
    if((vPosition.x() >= 0.) &&
       (vPosition.x() < vInterplanarPeriod)){
        return G4ThreeVector(fVectorEC->Value(vPosition.x()),0.,0.);
    }
    else{
        G4double vPositionX = vPosition.x() - fmod(vPosition.x(),vInterplanarPeriod) * vInterplanarPeriod;
        return G4ThreeVector(fVectorEC->Value(vPositionX),0.,0.);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector XVCrystalPlanarAnalytical::ComputePositionInUnitCell(G4ThreeVector vPosition, XPhysicalLattice* vLattice){
    G4VPhysicalVolume* vVolume = vLattice->GetVolume();

    G4double vInterplanarPeriod = GetXPhysicalLattice(vVolume)->ComputeInterplanarPeriod();
    
    G4double vPositionX = vPosition.x();
    
    if((vPositionX >= 0.) &&
       (vPositionX < vInterplanarPeriod)){
        return G4ThreeVector(vPositionX/vInterplanarPeriod,0.,0.);
    }
    else{
        vPositionX -= fmod(vPosition.x(),vInterplanarPeriod) * vInterplanarPeriod;
        return G4ThreeVector(vPositionX/vInterplanarPeriod,0.,0.);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XVCrystalPlanarAnalytical::ComputeMaximum(XPhysicalLattice* vLattice){
    unsigned int vPrecision = 1024;
    G4VPhysicalVolume* vVolume = vLattice->GetVolume();
    G4double vStep = GetXPhysicalLattice(vVolume)->ComputeInterplanarPeriod() / vPrecision;
    
    G4double vMaximum = -DBL_MAX;
    G4double vValue;
    
    for(unsigned int i=0;i<vPrecision;i++){
        if( (vValue = GetEC(G4ThreeVector(vStep * i,0.,0.),vLattice).x() ) > vMaximum) vMaximum = vValue;
    }
    return vMaximum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XVCrystalPlanarAnalytical::ComputeMinimum(XPhysicalLattice* vLattice){
    unsigned int vPrecision = 1024;
    G4VPhysicalVolume* vVolume = vLattice->GetVolume();
    G4double vStep = GetXPhysicalLattice(vVolume)->ComputeInterplanarPeriod() / vPrecision;
    
    G4double vMinimum = +DBL_MAX;
    G4double vValue;
    
    for(unsigned int i=0;i<vPrecision;i++){
        if( (vValue = GetEC(G4ThreeVector(vStep * i,0.,0.),vLattice).x() ) < vMinimum) vMinimum = vValue;
    }
    return vMinimum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XVCrystalPlanarAnalytical::PrintOnFile(const G4String& filename,XPhysicalLattice* vLattice,G4double vUnit){
    std::ofstream vFileOut;
    vFileOut.open(filename);
    vFileOut << "pos,val" << std::endl;
    
    G4int imax = 8192;
    G4double vXposition = 0.;
    G4double vInterplanarPeriod = vLattice->ComputeInterplanarPeriod();
    
    for(G4int i = 0;i<imax;i++){
        vXposition = double(i) / double(imax) * vInterplanarPeriod;
        vFileOut << vXposition / CLHEP::angstrom << "," << GetEC(G4ThreeVector(vXposition,0.,0.),vLattice).x() / vUnit << std::endl;
    }
    
    vFileOut.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XVCrystalPlanarAnalytical::InitializeVector(){
    G4int imax = 4096;
    G4double vXposition = 0.;
    G4double vInterplanarPeriod = fPhysicalLattice->ComputeInterplanarPeriod();
    
    fVectorEC = new G4PhysicsLinearVector(0,vInterplanarPeriod*(imax-1)/imax,imax);
    for(G4int i = 0;i<imax;i++){
        vXposition = double(i) / double(imax) * vInterplanarPeriod;
        fVectorEC->PutValue(i,ComputeEC(G4ThreeVector(vXposition,0.,0.),fPhysicalLattice).x());
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

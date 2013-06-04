//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//

#include "XCrystalIntegratedDensity.hh"

XCrystalIntegratedDensity::XCrystalIntegratedDensity(XVCrystalCharacteristic* vPotential,XVCrystalCharacteristic* vNucleiDensity,XVCrystalCharacteristic* vElectronDensity){
    
    SetNucleiDensity(vNucleiDensity);
    SetElectronDensity(vElectronDensity);
    SetPotential(vPotential);
    
    fTable.clear();
    
    fNumberOfPoints[0] = 256;
    fNumberOfPoints[1] = 1;
    fNumberOfPoints[2] = 1;
    
    fIntegrationPoints[0] = 1024;
    fIntegrationPoints[1] = 1;
    fIntegrationPoints[2] = 1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XCrystalIntegratedDensity::~XCrystalIntegratedDensity(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XCrystalIntegratedDensity::SetIntegrationPoints(unsigned int vIntegrationPoints,unsigned int vIndex){
    if(vIndex<3) {
        if(vIntegrationPoints > 0){
            fIntegrationPoints[vIndex] = vIntegrationPoints;
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

unsigned int XCrystalIntegratedDensity::GetIntegrationPoints(unsigned int vIndex){
    if(vIndex<3) {
        return fIntegrationPoints[vIndex];
    }
    else{
        return 0;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

unsigned int XCrystalIntegratedDensity::GetIntegrationPoints(){
    return fIntegrationPoints[0]*fIntegrationPoints[1]*fIntegrationPoints[2];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XCrystalIntegratedDensity::SetNumberOfPoints(unsigned int vNumberOfPoints,unsigned int vIndex){
    if(vIndex<3) {
        if(vNumberOfPoints > 0){
            fNumberOfPoints[vIndex] = vNumberOfPoints;
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

unsigned int XCrystalIntegratedDensity::GetNumberOfPoints(unsigned int vIndex){
    if(vIndex<3) {
        return fNumberOfPoints[vIndex];
    }
    else{
        return 0;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

unsigned int XCrystalIntegratedDensity::GetNumberOfPoints(){
    return fNumberOfPoints[0]*fNumberOfPoints[1]*fNumberOfPoints[2];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XCrystalIntegratedDensity::SetPotential(XVCrystalCharacteristic* vPotential){
    fPotential = vPotential;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XVCrystalCharacteristic* XCrystalIntegratedDensity::GetPotential(){
    return fPotential;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XCrystalIntegratedDensity::SetNucleiDensity(XVCrystalCharacteristic* vNucleiDensity){
    fNucleiDensity = vNucleiDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XVCrystalCharacteristic* XCrystalIntegratedDensity::GetNucleiDensity(){
    return fNucleiDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XCrystalIntegratedDensity::SetElectronDensity(XVCrystalCharacteristic* vElectronDensity){
    fElectronDensity = vElectronDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XVCrystalCharacteristic*  XCrystalIntegratedDensity::GetElectronDensity(){
    return fElectronDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XCrystalIntegratedDensity::SetXPhysicalLattice(XPhysicalLattice* vLattice){
    fLattice = vLattice;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XPhysicalLattice*  XCrystalIntegratedDensity::GetXPhysicalLattice(){
    return fLattice;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XCrystalIntegratedDensity::InitializeTable(){
    fPotentialMinimum = fPotential->GetMinimum(fLattice).x();
    fPotentialRange = fPotential->GetMaximum(fLattice).x() - fPotentialMinimum;
    
    for(unsigned int i=0;i<GetNumberOfPoints();i++){
        fTable.push_back(ComputeValue(fPotentialMinimum + fPotentialRange * G4double(i+1) / G4double(fNumberOfPoints[0])));
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XCrystalIntegratedDensity::GetValue(G4double vPotential){
    if(!HasBeenInitialized()) return -1.;
    else if(vPotential >= (fPotentialMinimum+fPotentialRange)) return 1.;
    else if(vPotential <= fPotentialMinimum) return 0.;
    else{
        unsigned int vIndex = int(((vPotential - fPotentialMinimum) / fPotentialRange ) * fNumberOfPoints[0]);
        vPotential -= fPotentialMinimum;
        vPotential /= fPotentialRange;
        return fTable.at(vIndex);
    }
    return 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool XCrystalIntegratedDensity::HasBeenInitialized(){
    if(fTable.size()==0) return false;
    else return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XCrystalIntegratedDensity::ComputeValue(G4double vPotentialInitial){
    
    unsigned int i1 = 0;
    
    G4ThreeVector vPositionTemp = G4ThreeVector(0.,0.,0.);
    G4double vDensity = 0.;
    
    G4double vInterplanarPeriod = fLattice->ComputeInterplanarPeriod();
    while(i1<fIntegrationPoints[0]){
        vPositionTemp.setX(G4double(G4double(i1)/G4double(fIntegrationPoints[0])*vInterplanarPeriod));
        if(fPotential->ComputeValue(vPositionTemp,fLattice).x() < vPotentialInitial){
            vDensity += fElectronDensity->ComputeValue(vPositionTemp,fLattice).x();
            vDensity += fNucleiDensity->ComputeValue(vPositionTemp,fLattice).x();
        }
        i1++;
    };
    vDensity *= vInterplanarPeriod;
    vDensity /= fIntegrationPoints[0];
    vDensity /= 2.;
    
    return vDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

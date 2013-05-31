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

XCrystalIntegratedDensity::XCrystalIntegratedDensity(XVCrystalCharacteristic* vNucleiDensity,XVCrystalCharacteristic* vElectronDensity){
    SetNucleiDensity(vNucleiDensity);
    SetElectronDensity(vElectronDensity);
    fTable.clear();
    fNumberOfTerms = 1024;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XCrystalIntegratedDensity::~XCrystalIntegratedDensity(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XCrystalIntegratedDensity::SetNumberOfPoints(unsigned int vNumberOfTerms){
    fNumberOfTerms = vNucleiDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

unsigned int XCrystalIntegratedDensity::GetNumberOfPoints(){
    return fNumberOfTerms;
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

void XCrystalIntegratedDensity::InitializeTable(XPhysicalLattice* vLattice){
vLattice
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XCrystalIntegratedDensity::FindIndex(G4ThreeVector vPosition) {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XCrystalIntegratedDensity::GetValue(G4ThreeVector vPosition){
    if(fTable.size() == 0) return -1.;
    else return fTable.at(FindIndex(vPosition));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

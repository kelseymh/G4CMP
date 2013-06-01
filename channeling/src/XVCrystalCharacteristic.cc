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

#include "XVCrystalCharacteristic.hh"

XVCrystalCharacteristic::XVCrystalCharacteristic(){
    fLatticeManager = XLatticeManager3::GetXLatticeManager();
    fThermalVibrationAmplitude = 0.075 * angstrom;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XVCrystalCharacteristic::~XVCrystalCharacteristic(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XPhysicalLattice* XVCrystalCharacteristic::GetXPhysicalLattice(G4VPhysicalVolume* vVolume)
{
    return fLatticeManager->GetXPhysicalLattice(vVolume);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XUnitCell* XVCrystalCharacteristic::GetXUnitCell(G4VPhysicalVolume* vVolume)
{
    return GetXPhysicalLattice(vVolume)->GetXUnitCell();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XLogicalLattice* XVCrystalCharacteristic::GetLogicalLattice(G4VPhysicalVolume* vVolume)
{
    return GetXPhysicalLattice(vVolume)->GetLogicalLattice();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XVCrystalCharacteristic::SetThermalVibrationAmplitude(G4double vThermalVibrationAmplitude){
    fThermalVibrationAmplitude = vThermalVibrationAmplitude;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XVCrystalCharacteristic::GetThermalVibrationAmplitude(){
    return fThermalVibrationAmplitude;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XVCrystalCharacteristic::ComputeTFScreeningRadius(XPhysicalLattice* vLattice){
    
    G4double vTFSR = Bohr_radius * 0.88534;

    vTFSR /= (std::pow(vLattice->GetXUnitCell()->GetBase(0)->GetElement()->GetZ(),0.333333333));

    //if(vLattice->GetParticleDefinition()->GetParticleName() == "proton"){
    //    vTFSR /= (std::pow(vLattice->GetMaterial()->GetZ(),0.333333333));
    //}
    //else{
    //    vTFSR /= (std::pow(vLattice->GetMaterial()->GetZ(),0.23) + std::pow(vLattice->GetParticleDefinition()->GetPDGCharge(),0.23));
    //}
    
    return vTFSR;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double ComputePositionInUnitCell{
    return -1.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double GetMaximum(){
    return (DBL_MAX);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double GetMinimum(){
    return (-DBL_MAX);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

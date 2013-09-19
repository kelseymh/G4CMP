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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "G4RunManager.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UniformMagField.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "A01DriftChamber.hh"
#include "A01DriftChamberHit.hh"

#include "XLatticeManager3.hh"

#include "XLogicalAtomicLattice.hh"
#include "XLogicalAtomicLatticeDiamond.hh"
#include "XLogicalBase.hh"
#include "XUnitCell.hh"

#include "XCrystalPlanarMolierePotential.hh"
#include "XCrystalPlanarMoliereElectricField.hh"
#include "XCrystalPlanarNucleiDensity.hh"
#include "XCrystalPlanarMoliereElectronDensity.hh"


#include "XCrystalPlanarMoliereTempPotential.hh"
#include "XCrystalPlanarMoliereElectricField.hh"
#include "XCrystalPlanarNucleiDensity.hh"
#include "XCrystalPlanarMoliereElectronDensity.hh"
#include "XCrystalCharacteristicArray.hh"
#include "XCrystalIntegratedDensityPlanar.hh"
#include "XCrystalIntegratedDensityHub.hh"
#include "XCrystalIntegratedDensityHub.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction():fWorldLogic(0),fXtalLogic(0){
    
    DefineMaterials();
    
    fWorldSizeY = 10. * CLHEP::meter;
    fWorldSizeXZ = 2. * CLHEP::meter;
    
    bSSD1 = false;
    bSSD2 = false;
    bSSD3 = false;
    fSSDSize = G4ThreeVector(1.92 * CLHEP::centimeter,0.06 * CLHEP::centimeter,1.92 * CLHEP::centimeter); //
    fSSD0XtalDistance = - 150 * CLHEP::centimeter;
    fSSD1XtalDistance = - 50 * CLHEP::centimeter;
    fSSD2XtalDistance = + 100 * CLHEP::centimeter;
    
    bSCI = false;
    fSCISize = G4ThreeVector(10. * CLHEP::centimeter,1. * CLHEP::centimeter,10. * CLHEP::centimeter); //
    fSCIRelativeDistance = 1. * CLHEP::centimeter;
    fSCIXtalDistance = 30. * CLHEP::centimeter;
    
    bRBS = false;
    fRBSSizeZ = 0.1 * CLHEP::centimeter;
    fRBSSizeR = 1 * CLHEP::centimeter;
    fRBSDistanceR = 10. * CLHEP::centimeter;
    fRBSAnglePhi = -135. * CLHEP::degree;
    fRBSAngleTheta = 0.;
    
    bXtal = false;
    fXtalAngle = G4ThreeVector(0.,0.,0.E-6 * CLHEP::radian);
    fXtalSize = G4ThreeVector(60. * CLHEP::millimeter, 2. * CLHEP::millimeter, 60. * CLHEP::millimeter);
    fXtalCellSize = G4ThreeVector(5.431 * CLHEP::angstrom, 5.431 * CLHEP::angstrom, 5.431 * CLHEP::angstrom);
    fXtalSize = G4ThreeVector(0. * CLHEP::millimeter, 0. * CLHEP::millimeter, 0. * CLHEP::millimeter);
    fXtalThermalVibrationAmplitude = 0.075 * CLHEP::angstrom;
    
    SetXtalMaterial("G4_Si");
    
    fMessenger = new DetectorConstructionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

DetectorConstruction::~DetectorConstruction(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void DetectorConstruction::DefineMaterials(){
    G4double const Torr = CLHEP::atmosphere/760.; // 1 Torr
    G4double z = 1.;
    G4double a = 1.01*CLHEP::g/CLHEP::mole;
    G4double density     = CLHEP::universe_mean_density;
    G4double pressure    = 1.E-5 * Torr;
    G4double temperature = 300.*CLHEP::kelvin;
    G4Material* Vacuum = new G4Material("Vacuum", z , a , density, kStateGas,temperature,pressure);
    
    fWorldMaterial = Vacuum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct(){
    
    ConstructWorld();
    
    if(bSCI) ConstructScintillators();
    
    if(bSSD1 || bSSD2 || bSSD3) ConstructSiliconStripDetectors();
    
    if(bRBS) ConstructRBSDetector();
    
    if(bXtal) ConstructXtalTarget();
    
    return fWorldPhysical;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructWorld(){
    
    if(fabs(fSSD2XtalDistance) > (fWorldSizeY/2.)){
        fWorldSizeY = fabs(fSSD2XtalDistance) * 3.;
    }
    if(fabs(fSSD1XtalDistance) > (fWorldSizeY/2.)){
        fWorldSizeY = fabs(fSSD1XtalDistance) * 3.;
    }
    if(fabs(fSSD0XtalDistance) > (fWorldSizeY/2.)){
        fWorldSizeY = fabs(fSSD0XtalDistance) * 3.;
    }
    
    fWorldSolid = new G4Box("World",fWorldSizeXZ/2.,fWorldSizeY/2.,fWorldSizeXZ/2.);
    
    fWorldLogic = new G4LogicalVolume(fWorldSolid,fWorldMaterial,"World");
    
    fWorldPhysical = new G4PVPlacement(0,G4ThreeVector(),fWorldLogic,"World",0,false,0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::AddSiliconStripDetectors(G4String vDetectorString) {
    if(vDetectorString == "0"){
        bSSD1 = true;
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
    }
    else if(vDetectorString == "1"){
        bSSD2 = true;
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
    }
    else if(vDetectorString == "2"){
        bSSD3 = true;
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
    }
    else if(vDetectorString == "all"){
        bSSD1 = true;
        bSSD2 = true;
        bSSD3 = true;
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructRBSDetector(){
    G4Material* Si = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");
    
    fRBSSolid = new G4Tubs("DetRBS", 0.,fRBSSizeR,fRBSSizeZ/2.,0,2*M_PI*CLHEP::radian);
    fRBSLogic = new G4LogicalVolume(fRBSSolid,Si,"DetRBS");
    
    
    G4ThreeVector vRBSDistance = G4ThreeVector(1.*CLHEP::centimeter,0.,0.);
    vRBSDistance.setR(fRBSDistanceR);
    vRBSDistance.rotate(G4ThreeVector(0,1,0), fRBSAngleTheta).rotate(G4ThreeVector(0,0,1), fRBSAnglePhi);
    G4RotationMatrix* vRotationMatrix = new G4RotationMatrix; // Rotates X and Z axes only
    vRotationMatrix->rotateX(M_PI/2 - fRBSAngleTheta);
    vRotationMatrix->rotateY(M_PI/2 + fRBSAnglePhi);
    fRBSPhysical = new G4PVPlacement(vRotationMatrix,vRBSDistance,fRBSLogic,"DetRBS",fWorldLogic,false,0);
    
    //    int vNumberOfDetectors = 32;
    //
    //    G4double vUnitAngle = M_PI * 3. / 2. / vNumberOfDetectors;
    //    for(int j1=0;j1<vNumberOfDetectors;j1++){
    //        G4ThreeVector vRBSDistance = G4ThreeVector(1.*CLHEP::centimeter,0.,0.);
    //        vRBSDistance.setR(fRBSDistanceR);
    //        vRBSDistance.rotate(G4ThreeVector(0,1,0), fRBSAngleTheta).rotate(G4ThreeVector(0,0,1), - M_PI/4 + j1 * vUnitAngle);
    //        G4RotationMatrix* vRotationMatrix = new G4RotationMatrix; // Rotates X and Z axes only
    //        vRotationMatrix->rotateZ( - M_PI/4 - j1 * vUnitAngle);
    //        vRotationMatrix->rotateX(M_PI/2 - fRBSAngleTheta);
    //        fRBSPhysical = new G4PVPlacement(vRotationMatrix,vRBSDistance,fRBSLogic,"DetRBS",fWorldLogic,false,j1);
    //    }
    
    G4String SDname;
    G4VSensitiveDetector* detectorRBS = new A01DriftChamber(SDname="/detectorRBS");
    G4SDManager::GetSDMpointer()->AddNewDetector(detectorRBS);
    fRBSLogic->SetSensitiveDetector(detectorRBS);
    
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructScintillators(){
    G4Material* NaI = G4NistManager::Instance()->FindOrBuildMaterial("G4_SODIUM_IODIDE");
    
    fSCISolid = new G4Box("NaISCI", fSCISize.x()/2.,fSCISize.y()/2.,fSCISize.z()/2.);
    fSCILogic = new G4LogicalVolume(fSCISolid,NaI,"NaISCI");
    
    G4double vDetDistX = (fSCIRelativeDistance / 2. + fSCISize.x() / 2.);
    
    new G4PVPlacement(0,G4ThreeVector(-vDetDistX,fSCIXtalDistance,0.),fSCILogic,"NaISCI",fWorldLogic,false,0);
    
    new G4PVPlacement(0,G4ThreeVector(+vDetDistX,fSCIXtalDistance,0.),fSCILogic,"NaISCI",fWorldLogic,false,1);
    
    G4String SDname;
    G4VSensitiveDetector* telescope = new A01DriftChamber(SDname="/scintillator");
    G4SDManager::GetSDMpointer()->AddNewDetector(telescope);
    fSCILogic->SetSensitiveDetector(telescope);
    
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSiliconStripDetectors(){
    G4Material* Si = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");
    
    fSSDSolid = new G4Box("SiSD", fSSDSize.x()/2.,fSSDSize.y()/2.,fSSDSize.z()/2.);
    fSSDLogic = new G4LogicalVolume(fSSDSolid,Si,"SiSD");
    
    if(bSSD1){
        new G4PVPlacement(0,G4ThreeVector(0.*CLHEP::centimeter,fSSD0XtalDistance,0.*CLHEP::centimeter),fSSDLogic,"SiSD",fWorldLogic,false,0);
    }
    if(bSSD2){
        new G4PVPlacement(0,G4ThreeVector(0.,fSSD1XtalDistance,0.),fSSDLogic,"SiSD",fWorldLogic,false,1);
    }
    if(bSSD3) {
        new G4PVPlacement(0,G4ThreeVector(0.*CLHEP::centimeter,fSSD2XtalDistance,0.*CLHEP::centimeter),fSSDLogic,"SiSD",fWorldLogic,false,2);
    }
    
    G4String SDname;
    G4VSensitiveDetector* telescope = new A01DriftChamber(SDname="/telescope");
    G4SDManager::GetSDMpointer()->AddNewDetector(telescope);
    fSSDLogic->SetSensitiveDetector(telescope);
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructXtalTarget(){
//    if(fXtalCurvatureRadius.x() != 0.){
//        fXtalSolid = new G4Tubs("Target",
//                                fXtalCurvatureRadius.x() - fXtalSize.x()/2,
//                                fXtalCurvatureRadius.x() + fXtalSize.x()/2,
//                                fXtalSize.z()/2.,
//                                fXtalAngle.z(),
//                                fXtalAngle.z() + fXtalSize.y()/fXtalCurvatureRadius.x());
//    }
//    else{
//        fXtalSolid = new G4Box("Target", fXtalSize.x()/2.,fXtalSize.y()/2.,fXtalSize.z()/2.);
//    }
    fXtalSolid = new G4Box("Target", fXtalSize.x()/2.,fXtalSize.y()/2.,fXtalSize.z()/2.);

    
    fXtalLogic = new G4LogicalVolume(fXtalSolid,fXtalMaterial,"Target");
    
    G4RotationMatrix* vRotationMatrix = new G4RotationMatrix; // Rotates X and Z axes only
    vRotationMatrix->rotateX(fXtalAngle.x());
    vRotationMatrix->rotateY(fXtalAngle.y());
    vRotationMatrix->rotateZ(fXtalAngle.z());
    
//    if(fXtalCurvatureRadius.x() != 0.){
//        fXtalPhysical = new G4PVPlacement(vRotationMatrix,G4ThreeVector(-fXtalCurvatureRadius.x(),0.,0.),fXtalLogic,"Target",fWorldLogic,false,0);
//    }
//    else{
//        fXtalPhysical = new G4PVPlacement(vRotationMatrix,G4ThreeVector(0.,0.,0.),fXtalLogic,"Target",fWorldLogic,false,0);
//    }
    
    fXtalPhysical = new G4PVPlacement(vRotationMatrix,G4ThreeVector(0.,0.,0.),fXtalLogic,"Target",fWorldLogic,false,0);

    //----------------------------------------
    // Create XLogicalLattice
    //----------------------------------------
    XLogicalLattice* logicalLattice = new XLogicalLattice();
    logicalLattice->SetScatteringConstant(3.67e-41*CLHEP::second*CLHEP::second*CLHEP::second);
    
    //----------------------------------------
    // Create XLogicalBase
    //----------------------------------------
    XLogicalAtomicLatticeDiamond *diamond_lattice = new XLogicalAtomicLatticeDiamond();
    G4Element* element = G4NistManager::Instance()->FindOrBuildElement(fXtalMaterial->GetZ());
    XLogicalBase *base = new XLogicalBase(element,diamond_lattice);
    
    //----------------------------------------
    // Create XUnitCell
    //----------------------------------------
    XUnitCell* myCell = new XUnitCell();
    myCell->SetSize(fXtalCellSize);
    myCell->AddBase(base);
    
    //----------------------------------------
    // Create XPhysicalLattice
    //----------------------------------------
    XPhysicalLattice* physicalLattice = new XPhysicalLattice(fXtalPhysical, logicalLattice);
    physicalLattice->SetUnitCell(myCell);
    physicalLattice->SetMillerOrientation(2,2,0);
    physicalLattice->SetLatticeOrientation(fXtalAngle.x(),fXtalAngle.y(),fXtalAngle.z());
    physicalLattice->SetThermalVibrationAmplitude(fXtalThermalVibrationAmplitude);
    physicalLattice->SetCurvatureRadius(fXtalCurvatureRadius);
    
    //----------------------------------------
    // Register XPhysicalLattice
    //----------------------------------------
    if(XLatticeManager3::GetXLatticeManager()->GetXPhysicalLattice(fXtalPhysical) != physicalLattice){
        XLatticeManager3::GetXLatticeManager()->RegisterLattice(physicalLattice);
    }
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldMaterial(const G4String& name){
    G4Material* vMaterial = G4Material::GetMaterial(name, false);
    
    if(!vMaterial){
        vMaterial = G4NistManager::Instance()->FindOrBuildMaterial(name);
    }
    
    if (vMaterial && vMaterial != fWorldMaterial) {
        G4cout << "DetectorConstructor::SetXtalMaterial() - New Xtal Material: " << vMaterial->GetName() << G4endl;
        fWorldMaterial = vMaterial;
        if(fWorldLogic){
            fWorldLogic->SetMaterial(vMaterial);
            G4RunManager::GetRunManager()->PhysicsHasBeenModified();
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String DetectorConstruction::GetWorldMaterial(){
    if(fXtalMaterial) {
        return fWorldMaterial->GetName();
    }
    return "";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.

void DetectorConstruction::SetXtalMaterial(const G4String& name){
    G4Material* vMaterial = G4Material::GetMaterial(name, false);
    
    if(!vMaterial){
        vMaterial = G4NistManager::Instance()->FindOrBuildMaterial(name);
    }
    
    if (vMaterial && vMaterial != fXtalMaterial) {
        G4cout << "DetectorConstructor::SetXtalMaterial() - New Xtal Material: " << vMaterial->GetName() << G4endl;
        fXtalMaterial = vMaterial;
        if(fXtalLogic){
            fXtalLogic->SetMaterial(vMaterial);
            G4RunManager::GetRunManager()->PhysicsHasBeenModified();
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String DetectorConstruction::GetXtalMaterial(){
    if(fXtalMaterial) {
        return fXtalMaterial->GetName();
    }
    return "";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSSDSize(G4ThreeVector size){
    if(fSSDSize != size) {
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
        fSSDSize = size;
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSSD0XtalDistance(G4double distance) {
    if(fSSD0XtalDistance != distance) {
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
        fSSD0XtalDistance = distance;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSSD1XtalDistance(G4double distance) {
    if(fSSD1XtalDistance != distance) {
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
        fSSD1XtalDistance = distance;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSSD2XtalDistance(G4double distance) {
    if(fSSD2XtalDistance != distance) {
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
        fSSD2XtalDistance = distance;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSCISize(G4ThreeVector size){
    if(fSCISize != size) {
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
        fSCISize = size;
    }
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSCIRelativeDistance(G4double distance) {
    if(fSCIRelativeDistance != distance) {
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
        fSCIRelativeDistance = distance;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSCIXtalDistance(G4double distance) {
    if(fSCIXtalDistance != distance) {
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
        fSCIXtalDistance = distance;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetRBSDistanceR(G4double distance) {
    if(fRBSDistanceR != distance) {
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
        fRBSDistanceR = distance;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetRBSAngleTheta(G4double angle) {
    if(fRBSAngleTheta != angle) {
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
        fRBSAngleTheta = angle;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetRBSAnglePhi(G4double angle) {
    if(fRBSAnglePhi != angle) {
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
        fRBSAnglePhi = angle;
        G4cout << "DetectorConstruction::SetRBSAnglePhi() - New Phi = " << fRBSAnglePhi/CLHEP::degree << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetXtalCurvatureRadius(G4ThreeVector cr) {
    if(fXtalCurvatureRadius != cr) {
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
        fXtalCurvatureRadius = cr;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetXtalSize(G4ThreeVector size) {
    if(fXtalSize != size) {
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
        fXtalSize = size;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetXtalAngle(G4ThreeVector angle) {
    if(fXtalAngle != angle) {
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
        fXtalAngle = angle;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetXtalCellSize(G4ThreeVector cellsize) {
    if(fXtalCellSize != cellsize) {
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
        fXtalCellSize = cellsize;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetXtalCellAngle(G4ThreeVector cellangle) {
    if(fXtalCellAngle != cellangle) {
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
        fXtalCellAngle = cellangle;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetXtalThermalVibrationAmplitude(G4double thermvibr) {
    if(fXtalThermalVibrationAmplitude != thermvibr) {
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
        fXtalThermalVibrationAmplitude = thermvibr;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo..

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
//:fMessenger(0),
//bXtal(0),fXtalCurvatureRadius(0),fXtalMaterial(0),fXtalAngle(0),fXtalSize(0),
//fXtalCellSize(0),fXtalCellAngle(0),fXtalSolid(0),fXtalLogic(0),fXtalPhysical(0),
//bRBS(0),fRBSDistanceR(0),fRBSAngleTheta(0),fRBSAnglePhi(0),
//fRBSSizeZ(0),fRBSSizeR(0),fRBSSolid(0),fRBSLogic(0),fRBSPhysical(0),
//bSCI(0),fSCIXtalDistance(0),fSCIRelativeDistance(0),
//fSCISizeXZ(0),fSCISizeY(0),fSCISolid(0),fSCILogic(0),fSCIPhysical(0),
//bSSD(0),fSSD0XtalDistance(0),fSSD1XtalDistance(0),fSSD2XtalDistance(0),
//fSSDSizeXZ(0),fSSDSizeY(0),fSSDSolid(0),fSSDLogic(0),fSSDPhysical(0),
//fWorldSizeXZ(0),fWorldSizeY(0),fWorldSolid(0),fWorldLogic(0),fWorldPhysical(0)
{

    DefineMaterials();
    
    fWorldSizeY = 0.5 * m;
    fWorldSizeXZ = 0.5 * m;
    
    bSSD = false;
    fSSDSizeXZ = 1.92 * 10. * cm; // originally 1.92 cm
    fSSDSizeY = 0.06 * cm;

    fSSD0XtalDistance = - 15 * cm;
    fSSD1XtalDistance = - 5 * cm;
    fSSD2XtalDistance = + 10 * cm;

    bSCI = false;
    fSCIRelativeDistance = 1. * cm;
    fSCIXtalDistance = 30. * cm;
    fSCISizeXZ = 10. * cm;
    fSCISizeY = 1. * cm;
    
    bRBS = false;
    fRBSSizeZ = 0.1 * cm;
    fRBSSizeR = 1 * cm;
    fRBSDistanceR = 10. * cm;
    fRBSAnglePhi = -135. * degree;
    fRBSAngleTheta = 0.;

    bXtal = false;
    fXtalAngle = G4ThreeVector(0.,0.,0.E-6 * radian);
    fXtalSize = G4ThreeVector(60. * mm, 2. * mm, 60. * mm);
    fXtalCellSize = G4ThreeVector(5.431 * angstrom, 5.431 * angstrom, 5.431 * angstrom);
    fXtalSize = G4ThreeVector(0. * mm, 0. * mm, 0. * mm);
    
    SetXtalMaterial("G4_Si");
    
    fMessenger = new DetectorConstructionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

DetectorConstruction::~DetectorConstruction(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void DetectorConstruction::DefineMaterials(){
    G4Material* Si = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");
    G4Material* Ge = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ge");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct(){
    
    ConstructWorld();
        
    if(bSCI) ConstructScintillators();
    
    if(bSSD) ConstructSiliconStripDetectors();
    
    if(bRBS) ConstructRBSDetector();
    
    if(bXtal) ConstructXtalTarget();
    
    return fWorldPhysical;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructWorld(){
    G4Material* Vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    
    if(fabs(fSSD2XtalDistance) > fWorldSizeY/2.){
       fWorldSizeY = fabs(fSSD2XtalDistance) * 2.6;
    }
    if(fabs(fSSD1XtalDistance) > fWorldSizeY/2.){
        fWorldSizeY = fabs(fSSD1XtalDistance) * 2.6;
    }
    if(fabs(fSSD0XtalDistance) > fWorldSizeY/2.){
        fWorldSizeY = fabs(fSSD0XtalDistance) * 2.6;
    }
    
    fWorldSolid = new G4Box("World",fWorldSizeXZ/2.,fWorldSizeY/2.,fWorldSizeXZ/2.);

    fWorldLogic = new G4LogicalVolume(fWorldSolid,Vacuum,"World");
    
    fWorldPhysical = new G4PVPlacement(0,G4ThreeVector(),fWorldLogic,"World",0,false,0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructRBSDetector(){
    G4Material* Si = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");
    
    G4ThreeVector vRBSDistance = G4ThreeVector(1.*cm,0.,0.);
    
    vRBSDistance.setR(fRBSDistanceR);
    vRBSDistance.rotate(G4ThreeVector(0,1,0), fRBSAngleTheta).rotate(G4ThreeVector(0,0,1), fRBSAnglePhi);
        
    fRBSSolid = new G4Tubs("DetRBS", 0.,fRBSSizeR,fRBSSizeZ/2.,0,2*M_PI*radian);
    fRBSLogic = new G4LogicalVolume(fRBSSolid,Si,"DetRBS");
        
    G4RotationMatrix* vRotationMatrix = new G4RotationMatrix; // Rotates X and Z axes only
    vRotationMatrix->rotateX(M_PI/2 - fRBSAngleTheta);
    vRotationMatrix->rotateY(M_PI/2 + fRBSAnglePhi);

    fRBSPhysical = new G4PVPlacement(vRotationMatrix,vRBSDistance,fRBSLogic,"DetRBS",fWorldLogic,false,0);
    
    G4String SDname;
    G4VSensitiveDetector* detectorRBS = new A01DriftChamber(SDname="/detectorRBS");
    G4SDManager::GetSDMpointer()->AddNewDetector(detectorRBS);
    fRBSLogic->SetSensitiveDetector(detectorRBS);
    
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructScintillators(){
    G4Material* NaI = G4NistManager::Instance()->FindOrBuildMaterial("G4_SODIUM_IODIDE");
     
    fSCISolid = new G4Box("NaISCI", fSCISizeXZ/2.,fSCISizeY/2.,fSCISizeXZ/2.);
    fSCILogic = new G4LogicalVolume(fSCISolid,NaI,"NaISCI");
    
    G4double vDetDistX = (fSCIRelativeDistance / 2. + fSCISizeXZ / 2.);

    fSCIPhysical = new G4PVPlacement(0,G4ThreeVector(-vDetDistX,fSCIXtalDistance,0.),fSCILogic,"NaISCI",fWorldLogic,false,0);
    
    fSCIPhysical = new G4PVPlacement(0,G4ThreeVector(+vDetDistX,fSCIXtalDistance,0.),fSCILogic,"NaISCI",fWorldLogic,false,1);

    G4String SDname;
    G4VSensitiveDetector* telescope = new A01DriftChamber(SDname="/scintillator");
    G4SDManager::GetSDMpointer()->AddNewDetector(telescope);
    fSCILogic->SetSensitiveDetector(telescope);
    
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSiliconStripDetectors(){
    G4Material* Si = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");
    
    fSSDSolid = new G4Box("SiSD", fSSDSizeXZ/2.,fSSDSizeY/2.,fSSDSizeXZ/2.);
    fSSDLogic = new G4LogicalVolume(fSSDSolid,Si,"SiSD");
    
    fSSDPhysical = new G4PVPlacement(0,G4ThreeVector(0.*cm,fSSD0XtalDistance,0.*cm),fSSDLogic,"SiSD",fWorldLogic,false,0);
    
    fSSDPhysical = new G4PVPlacement(0,G4ThreeVector(0.*cm,fSSD1XtalDistance,0.*cm),fSSDLogic,"SiSD",fWorldLogic,false,1);
    
    fSSDPhysical = new G4PVPlacement(0,G4ThreeVector(0.*cm,fSSD2XtalDistance,0.*cm),fSSDLogic,"SiSD",fWorldLogic,false,2);

    G4String SDname;
    G4VSensitiveDetector* telescope = new A01DriftChamber(SDname="/telescope");
    G4SDManager::GetSDMpointer()->AddNewDetector(telescope);
    fSSDLogic->SetSensitiveDetector(telescope);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructXtalTarget(){
    fXtalSolid = new G4Box("Target", fXtalSize.x()/2.,fXtalSize.y()/2.,fXtalSize.z()/2.);
    fXtalLogic = new G4LogicalVolume(fXtalSolid,fXtalMaterial,"Target");

    G4RotationMatrix* vRotationMatrix = new G4RotationMatrix; // Rotates X and Z axes only
    vRotationMatrix->rotateX(fXtalAngle.x());
    vRotationMatrix->rotateY(fXtalAngle.y());
    vRotationMatrix->rotateZ(fXtalAngle.z());

    fXtalPhysical = new G4PVPlacement(vRotationMatrix,G4ThreeVector(0.,1.*cm,0.),fXtalLogic,"Target",fWorldLogic,false,0);
    
    //----------------------------------------
    // Create XLogicalLattice
    //----------------------------------------
    XLogicalLattice* logicalLattice = new XLogicalLattice();
    logicalLattice->SetScatteringConstant(3.67e-41*s*s*s);
    
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
    physicalLattice->SetLatticeOrientation(fXtalAngle.x(),fXtalAngle.z());
    physicalLattice->SetThermalVibrationAmplitude(0.075*angstrom);
    physicalLattice->SetCurvatureRadius(fXtalCurvatureRadius);
    
    //----------------------------------------
    // Register XPhysicalLattice
    //----------------------------------------
    if(XLatticeManager3::GetXLatticeManager()->GetXPhysicalLattice(fXtalPhysical) != physicalLattice){
        XLatticeManager3::GetXLatticeManager()->RegisterLattice(physicalLattice);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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
        G4cout << "DetectorConstruction::SetRBSAnglePhi() - New Phi = " << fRBSAnglePhi/deg << G4endl;
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

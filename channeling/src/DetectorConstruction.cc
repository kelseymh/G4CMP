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

DetectorConstruction::DetectorConstruction(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

DetectorConstruction::~DetectorConstruction(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct(){
    AddWorld();
    AddCrystalTarget();
    AddScintillators();
    AddSiliconStripDetectors();
    return fWorldPhysical;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::AddWorld(){
    G4Material* Vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    
    fWorldSizeY = 2.4 * m;
    fWorldSizeXZ = 1. * m;
    
    fWorldSolid = new G4Box("World", fWorldSizeXZ/2.,fWorldSizeY/2.,fWorldSizeXZ/2.);
    fWorldLogic = new G4LogicalVolume(fWorldSolid,Vacuum,"World");
    fWorldPhysical = new G4PVPlacement(0,G4ThreeVector(0,0,0),fWorldLogic,"World",0,false,0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::AddScintillators(){
    G4Material* Si = G4NistManager::Instance()->FindOrBuildMaterial("G4_SODIUM_IODIDE");
    
    fSCISizeXZ = 100. * mm;
    fSCISizeY = 0.001 * mm;
    
    G4double vSCI_CryDistance = 60. * cm;
    G4double vSCIXDistance = 10. * mm;

    fSCISolid = new G4Box("NaISCI", fSCISizeXZ/2.,fSCISizeY/2.,fSCISizeXZ/2.);
    fSCILogic = new G4LogicalVolume(fSCISolid,Si,"NaISCI");
    
    for(int j1=0;j1<2;j1++)
    {
        G4double vDetDistX = double(j1 - 1) * (vSCIXDistance + fSCISizeXZ/2.*2.);
        vDetDistX += (vSCIXDistance / 2. + fSCISizeXZ / 2.);
        fSCIPhysical = new G4PVPlacement(0,G4ThreeVector(vDetDistX,vSCI_CryDistance,0.),fSCILogic,"NaISCI",fWorldLogic,false,j1);
    }
    
    G4String SDname;
    G4VSensitiveDetector* telescope = new A01DriftChamber(SDname="/scintillator");
    G4SDManager::GetSDMpointer()->AddNewDetector(telescope);
    fSCILogic->SetSensitiveDetector(telescope);
    
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::AddSiliconStripDetectors(){
    G4Material* Si = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");
    
    fSSDSizeXZ = 1.92 * cm;
    fSSDSizeY = 0.6 * mm;
    
    G4double vSD_CryDistance = 1. * cm;
    G4double vSDXDistance = 70. * cm;

    fSSDSolid = new G4Box("SiSD", fSSDSizeXZ/2.,fSSDSizeY/2.,fSSDSizeXZ/2.);
    fSSDLogic = new G4LogicalVolume(fSSDSolid,Si,"SiSD");
    
    for(int j1=0;j1<3;j1++)
    {
        G4double vDetDistY = double(j1 - 1) * vSDXDistance + fSSDSizeY/2. - (vSD_CryDistance + fXtalSizeY/2.);;
        fSSDPhysical = new G4PVPlacement(0,G4ThreeVector(0.,vDetDistY,0.),fSSDLogic,"SiSD",fWorldLogic,false,j1);
    }
    
    G4String SDname;
    G4VSensitiveDetector* telescope = new A01DriftChamber(SDname="/telescope");
    G4SDManager::GetSDMpointer()->AddNewDetector(telescope);
    fSSDLogic->SetSensitiveDetector(telescope);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::AddCrystalTarget(){
    G4Material* Si = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");

    fXtalSizeY = 2. * mm;
    fXtalSizeXZ = 6. * cm;
    
    fXtalSolid = new G4Box("Target", fXtalSizeXZ/2.,fXtalSizeY/2.,fXtalSizeXZ/2.);
    fXtalLogic = new G4LogicalVolume(fXtalSolid, Si,"Target");
    fXtalPhysical = new G4PVPlacement(0,G4ThreeVector(),fXtalLogic,"Target",fWorldLogic,false,0);
    
    //----------------------------------------
    // Create XLogicalLattice
    //----------------------------------------
    XLogicalLattice* logicalLattice = new XLogicalLattice();
    logicalLattice->SetScatteringConstant(3.67e-41*s*s*s);
    
    //----------------------------------------
    // Create XLogicalBase
    //----------------------------------------
    XLogicalAtomicLatticeDiamond *diamond_lattice = new XLogicalAtomicLatticeDiamond();
    G4Element* element = G4NistManager::Instance()->FindOrBuildElement(14);
    XLogicalBase *base = new XLogicalBase(element,diamond_lattice);

    //----------------------------------------
    // Create XUnitCell
    //----------------------------------------
    XUnitCell* myCell = new XUnitCell();
    myCell->SetSize(G4ThreeVector(5.431 * angstrom, 5.431 * angstrom, 5.431 * angstrom));
    myCell->AddBase(base);
    
    //----------------------------------------
    // Create XPhysicalLattice
    //----------------------------------------
    XPhysicalLattice* physicalLattice = new XPhysicalLattice(fXtalPhysical, logicalLattice);
    physicalLattice->SetUnitCell(myCell);
    physicalLattice->SetMillerOrientation(2,2,0);
    physicalLattice->SetLatticeOrientation(0.,0.E-6 * radian);

    //----------------------------------------
    // Register XPhysicalLattice
    //----------------------------------------
    XLatticeManager3::GetXLatticeManager()->RegisterLattice(physicalLattice);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

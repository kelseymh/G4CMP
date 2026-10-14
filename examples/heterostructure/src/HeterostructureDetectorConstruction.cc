/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file HeterostructureDetectorConstruction.cc 
/// \brief Implementation of the HeterostructureDetectorConstruction class
/// Constructs heterostructure toy geometry with germanium on silicon.

// 20260815  Selby Q. Dang (Stanford/SLAC)

#include "HeterostructureDetectorConstruction.hh"
#include "HeterostructureConfigManager.hh"
#include "HeterostructureSensitivity.hh"
#include "G4CMPLogicalBorderSurface.hh"
#include "G4CMPLogicalSkinSurface.hh"
#include "G4CMPPhononElectrode.hh"
#include "G4CMPSurfaceProperty.hh"
#include "G4CMPUtils.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4GeometryManager.hh"
#include "G4LatticeLogical.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RunManager.hh"
#include "G4SolidStore.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

HeterostructureDetectorConstruction::HeterostructureDetectorConstruction()
  : fLiquidHelium(0), fGermanium(0), fSilicon(0), fWorldPhys(0), 
    fGeSiInterface(0), fSiGeInterface(0), fGeVacInterface(0), fSiVacInterface(0), 
    fConstructed(false) {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

HeterostructureDetectorConstruction::~HeterostructureDetectorConstruction() {
  delete fGeSiInterface;
  delete fSiGeInterface;
  delete fGeVacInterface;
  delete fSiVacInterface;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* HeterostructureDetectorConstruction::Construct() {
  if (fConstructed) {
    if (!G4RunManager::IfGeometryHasBeenDestroyed()) {
      // Run manager hasn't cleaned volume stores. This code shouldn't execute
      G4GeometryManager::GetInstance()->OpenGeometry();
      G4PhysicalVolumeStore::GetInstance()->Clean();
      G4LogicalVolumeStore::GetInstance()->Clean();
      G4SolidStore::GetInstance()->Clean();
    }
    // Have to completely remove all lattices to avoid warning on reconstruction
    G4LatticeManager::GetLatticeManager()->Reset();
    // Clear all LogicalSurfaces
    // NOTE: No need to redefine the G4CMPSurfaceProperties
    G4CMPLogicalBorderSurface::CleanSurfaceTable();
  }

  DefineMaterials();
  SetupGeometry();
  fConstructed = true;

  return fWorldPhys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HeterostructureDetectorConstruction::DefineMaterials()
{ 
  G4NistManager* nistManager = G4NistManager::Instance();

  fLiquidHelium = nistManager->FindOrBuildMaterial("G4_AIR"); // to be corrected
  fGermanium = nistManager->FindOrBuildMaterial("G4_Ge");
  fSilicon = nistManager->FindOrBuildMaterial("G4_Si");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HeterostructureDetectorConstruction::SetupGeometry()
{
  //     
  // World
  //
  G4VSolid* worldSolid = new G4Box("World",16.*cm,16.*cm,16.*cm);
  G4LogicalVolume* worldLogical =
    new G4LogicalVolume(worldSolid,fLiquidHelium,"World");
  worldLogical->SetUserLimits(new G4UserLimits(10*mm, DBL_MAX, DBL_MAX, 0, 0));
  fWorldPhys = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"World",0,
                                 false,0);
  
  //
  // Germanium layer -- should change to SiGe in the future!
  // 
  G4VSolid* fGermaniumSolid = new G4Box("fGermaniumSolid", 10.*mm, 10.*mm, 5.*mm);
  G4LogicalVolume* fGermaniumLogical =
    new G4LogicalVolume(fGermaniumSolid,fGermanium,"fGermaniumLogical");
  G4VPhysicalVolume* GePhys =
    new G4PVPlacement(0,G4ThreeVector(),fGermaniumLogical,"fGermaniumPhysical",
                      worldLogical,false,0);

  //
  //Germanium lattice information
  //
  // G4LatticeManager gives physics processes access to lattices by volume
  G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();
  G4LatticeLogical* GeLogical = LM->LoadLattice(fGermanium, "Ge");

  // G4LatticePhysical assigns G4LatticeLogical a physical orientation
  G4LatticePhysical* GePhysical = new G4LatticePhysical(GeLogical);
  GePhysical->SetMillerOrientation(1,0,0);
  LM->RegisterLattice(GePhys, GePhysical);

  //
  // Silicon layer
  //
  G4VSolid* fSiliconSolid = new G4Box("fSiliconSolid",10.*mm,10.*mm,5.*mm);
  G4LogicalVolume* fSiliconLogical =
    new G4LogicalVolume(fSiliconSolid,fSilicon,"fSiliconLogical");
  G4VPhysicalVolume* SiPhys = new G4PVPlacement(0,
    G4ThreeVector(0.,0.,10.*mm), fSiliconLogical, "fSiliconPhysical",
    worldLogical,false,0);

  //
  //Silicon lattice information
  //
  G4LatticeLogical* SiLogical = LM->LoadLattice(fSilicon, "Si");
  // G4LatticePhysical assigns G4LatticeLogical a physical orientation
  G4LatticePhysical* SiPhysical = new G4LatticePhysical(SiLogical);
  SiPhysical->SetMillerOrientation(1,0,0);
  LM->RegisterLattice(SiPhys, SiPhysical);

  double qAbsProbGeSi  = HeterostructureConfigManager::GetQAbsProbGeSi();
  double qReflProbGeSi = HeterostructureConfigManager::GetQReflProbGeSi();
  double qAbsProbSiGe  = HeterostructureConfigManager::GetQAbsProbSiGe();
  double qReflProbSiGe = HeterostructureConfigManager::GetQReflProbSiGe();


  if (!fConstructed) {

    fGeSiInterface = new G4CMPSurfaceProperty("fGeSiInterface", 
                                          qAbsProbGeSi, qReflProbGeSi, 0.0, 0.0,
                                          0.0, 0.0, 0.0, 0.0, 
                                          0.0, 0.0);

    fSiGeInterface = new G4CMPSurfaceProperty("fSiGeInterface", 
                                          qAbsProbSiGe, qReflProbSiGe, 0.0, 0.0,
					  	                            0.0, 0.0, 0.0, 0.0, 
                                          0.0, 0.0);

    fGeVacInterface = new G4CMPSurfaceProperty("GeVacInterface", 
                                            1.0, 1.0, 0.0, 0.0,
					    	                            1.0, 1.0, 0.0, 0.0,
                                            1.0, 1.0);

    fSiVacInterface = new G4CMPSurfaceProperty("SiVacInterface", 
                                            1.0, 1.0, 0.0, 0.0,
                                            1.0, 1.0, 0.0, 0.0,      
                                            1.0, 1.0);
  }

  //
  // Border surfaces
  //
  new G4CMPLogicalBorderSurface("GeSi", GePhys, SiPhys,
        fGeSiInterface);
  new G4CMPLogicalBorderSurface("SiGe", SiPhys, GePhys, 
        fSiGeInterface);

  // new G4CMPLogicalSkinSurface("SiLogicalSkinSurface", fSiliconLogical,
  //          fSiVacInterface);
  // new G4CMPLogicalSkinSurface("GeLogicalSkinSurface", fGermaniumLogical,
  //          fGeVacInterface);

  new G4CMPLogicalBorderSurface("GeVac", GePhys, fWorldPhys,
				fGeVacInterface);
  new G4CMPLogicalBorderSurface("VacGe", fWorldPhys, GePhys, 
        fGeVacInterface);

  new G4CMPLogicalBorderSurface("SiVac", SiPhys, fWorldPhys,
        fSiVacInterface);
  new G4CMPLogicalBorderSurface("VacSi", fWorldPhys, SiPhys, 
        fSiVacInterface);

  //                                        
  // Visualization attributes
  //
  worldLogical->SetVisAttributes(G4VisAttributes::GetInvisible());
  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  G4VisAttributes* siVisAtt= new G4VisAttributes(G4Colour(1.0,0,1.0));
  siVisAtt->SetVisibility(true);
  fGermaniumLogical->SetVisAttributes(simpleBoxVisAtt);
  fSiliconLogical->SetVisAttributes(siVisAtt);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HeterostructureDetectorConstruction::ConstructSDandField()
{
  
  G4SDManager* sdman = G4SDManager::GetSDMpointer();

  G4VSensitiveDetector* geSD = sdman->FindSensitiveDetector("GermaniumElectrode", false);
  if (!geSD){ // Only create detector if one doesn't exist.
    geSD = new HeterostructureSensitivity("GermaniumElectrode");
    sdman->AddNewDetector(geSD);
  }
  SetSensitiveDetector("fGermaniumLogical",geSD);

  G4VSensitiveDetector* siSD = sdman->FindSensitiveDetector("SiliconElectrode", false);
  if (!siSD) {
    siSD = new HeterostructureSensitivity("SiliconElectrode");
    sdman->AddNewDetector(siSD);
  }
  SetSensitiveDetector("fSiliconLogical", siSD);
}

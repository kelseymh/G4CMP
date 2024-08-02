/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/


//
//20240110 Israel Hernandez -- Illinois Institute of Technology


#include "Caustic_PhononDetectorConstruction.hh"
#include "Caustic_PhononSensitivity.hh"
#include "G4CMPLogicalBorderSurface.hh"
#include "G4CMPPhononElectrode.hh"
#include "G4CMPSurfaceProperty.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4GeometryManager.hh"
#include "G4LatticeLogical.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SolidStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4Tubs.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
#include "Caustic_PhononConfigManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Caustic_PhononDetectorConstruction::Caustic_PhononDetectorConstruction()
  : fLiquidHelium(0),fBolometer(0),fOxigen(0),fSubstrate(0),
    fWorldPhys(0), topSurfProp(0), wallSurfProp(0),
    electrodeSensitivity(0), fConstructed(false) {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Caustic_PhononDetectorConstruction::~Caustic_PhononDetectorConstruction() {
  delete topSurfProp;

  delete wallSurfProp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* Caustic_PhononDetectorConstruction::Construct()
{
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

  Caustic_DefineMaterials();
  Caustic_SetupGeometry();
  fConstructed = true;

  return fWorldPhys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Caustic_PhononDetectorConstruction::Caustic_DefineMaterials()
{
  G4NistManager* nistManager = G4NistManager::Instance();

  fLiquidHelium = nistManager->FindOrBuildMaterial("G4_AIR"); // to be corrected
  fBolometer = nistManager->FindOrBuildMaterial("G4_Al");

//This lines works for Sapphire
//
 fOxigen = nistManager->FindOrBuildMaterial("G4_O");
fSubstrate = new G4Material("fSubstrate", 3.98*g/cm3, 2);
fSubstrate->AddElement(nistManager->FindOrBuildElement("Al"), 2);
fSubstrate->AddElement(nistManager->FindOrBuildElement("O"), 3);
//
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Caustic_PhononDetectorConstruction::Caustic_SetupGeometry()
{
    bool checkOverlaps = true;
  //
  // World
  //
  G4VSolid* worldSolid = new G4Box("World",6.*cm,6.*cm,6.*cm);
  G4LogicalVolume* worldLogical =
    new G4LogicalVolume(worldSolid,fLiquidHelium,"World");
  worldLogical->SetUserLimits(new G4UserLimits(10*mm, DBL_MAX, DBL_MAX, 0, 0));
  fWorldPhys = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"World",0,
                                 false,0);

// Substrate
 G4VSolid* fSubstrateSolid= new G4Box("fSubstrateSolid",0.2*cm,0.2*cm,0.2*cm);

  G4LogicalVolume* fSubstrateLogical =
    new G4LogicalVolume(fSubstrateSolid,fSubstrate,"fSubstrateLogical");
  G4VPhysicalVolume* SubstratePhys =
    new G4PVPlacement(0,G4ThreeVector(),fSubstrateLogical,"fSubstratePhysical",
                      worldLogical,false,0,checkOverlaps);

  //
  //Substrate lattice information
  //

  // G4LatticeManager gives physics processes access to lattices by volume
  //Al2O3 (Sapphire ) is the name of the folder where all the physical properties of the  Substrate are saved
  G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();
  //If you want to see the phonon caustic of other materials only change the name of the latttice .
  //Example   G4LatticeLogical* GeLogical = LM->LoadLattice(fSubstrate, "Si"); // Here you will see the phonons caustic of Silicon
  G4LatticeLogical* SubstrateLogical = LM->LoadLattice(fSubstrate, "Al2O3");
  //G4LatticeLogical* SubstrateLogical = LM->LoadLattice(fSubstrate, "GaAs");// For GaAs

  G4LatticePhysical* SubstratePhysical = new G4LatticePhysical(SubstrateLogical);


  SubstratePhysical->SetMillerOrientation(0,0,0);


  LM->RegisterLattice(SubstratePhys, SubstratePhysical);
  //
  // Aluminum . This is where phonon hits are registered (Bolometer)
  //


    G4VSolid* fBolometerSolid= new G4Box("fBolometerSolid",0.2*cm,0.2*cm,0.001*cm);

  G4LogicalVolume* fBolometerLogical =
    new G4LogicalVolume(fBolometerSolid,fBolometer,"fBolometerLogical");
  G4VPhysicalVolume* fBolometerPhysical = new G4PVPlacement(0,
    G4ThreeVector(0.,0.,0.201*cm), fBolometerLogical, "fBolometerPhysical",
    worldLogical,false,0,checkOverlaps);

  //
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  if (!electrodeSensitivity)
    electrodeSensitivity = new Caustic_PhononSensitivity("PhononElectrode");
  SDman->AddNewDetector(electrodeSensitivity);
  fSubstrateLogical->SetSensitiveDetector(electrodeSensitivity);

  //
  // surface between bolometer and Substrate determines phonon reflection/absorption
  // For the Phonon Caustics we want total phonon absorption for the bolometer and the Substrate walls.
  //
  if (!fConstructed) {


    topSurfProp = new G4CMPSurfaceProperty("TopAlSurf", 1.0, 0.0, 0.0, 0.0,
					  	        1.0, 1.0, 0.0, 0.0);




    wallSurfProp = new G4CMPSurfaceProperty("WallSurf", 0.0, 1.0, 0.0, 0.0,
					    	          1.0, 1.0, 0.0, 0.0);


  }

  //
  // Separate surfaces for sensors vs. bare sidewall
  //
  new G4CMPLogicalBorderSurface("detTop", SubstratePhys, fBolometerPhysical,
				topSurfProp);

  new G4CMPLogicalBorderSurface("detWall", SubstratePhys, fWorldPhys,
				wallSurfProp);

  //
  // Visualization attributes
  //
  worldLogical->SetVisAttributes(G4VisAttributes::Invisible);
  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));// Gray Color
  G4VisAttributes* simpleDetectorAtt= new G4VisAttributes(G4Colour(0.0,0.0,1.0));// Blue Color
  simpleBoxVisAtt->SetVisibility(true);
  fSubstrateLogical->SetVisAttributes(simpleBoxVisAtt);
  fBolometerLogical->SetVisAttributes(simpleDetectorAtt);
}


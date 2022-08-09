/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file exoticphysics/phonon/src/PhononDetectorConstruction.cc \brief
/// Implementation of the PhononDetectorConstruction class
//
// $Id: a2016d29cc7d1e75482bfc623a533d20b60390da $
//
// 20140321  Drop passing placement transform to G4LatticePhysical
// 20211207  Replace G4Logical*Surface with G4CMP-specific versions.
// 20220809  [ For M. Hui ] -- Add frequency dependent surface properties.

#include "PhononDetectorConstruction.hh"
#include "PhononSensitivity.hh"
#include "G4CMPElectrodeSensitivity.hh"
#include "G4CMPLogicalBorderSurface.hh"
#include "G4CMPSurfaceProperty.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4FieldManager.hh"
#include "G4GeometryManager.hh"
#include "G4LatticeLogical.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4CMPLogicalBorderSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SolidStore.hh"
#include "G4Sphere.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4Tubs.hh"
#include "G4UniformMagField.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PhononDetectorConstruction::PhononDetectorConstruction()
  : fLiquidHelium(0), fGermanium(0), fAluminum(0), fTungsten(0),
    fWorldPhys(0), topSurfProp(0), botSurfProp(0), wallSurfProp(0),
    electrodeSensitivity(0), fConstructed(false), fIfField(true) {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PhononDetectorConstruction::~PhononDetectorConstruction() {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* PhononDetectorConstruction::Construct()
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

  DefineMaterials();
  SetupGeometry();
  fConstructed = true;

  return fWorldPhys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhononDetectorConstruction::DefineMaterials()
{ 
  G4NistManager* nistManager = G4NistManager::Instance();

  fLiquidHelium = nistManager->FindOrBuildMaterial("G4_AIR"); // to be corrected
  fGermanium = nistManager->FindOrBuildMaterial("G4_Ge");
  fAluminum = nistManager->FindOrBuildMaterial("G4_Al");
  fTungsten = nistManager->FindOrBuildMaterial("G4_W");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhononDetectorConstruction::SetupGeometry()
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
  // Germanium cylinder - this is the volume in which we will propagate phonons
  //  
  G4VSolid* fGermaniumSolid = new G4Tubs("fGermaniumSolid",0.*cm,3.81*cm,
                                         1.27*cm, 0.*deg, 360.*deg);
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

  // NOTE:  Above registration can also be done in single step:
  // G4LatticlePhysical* GePhysical = LM->LoadLattice(GePhys, "Ge");

  //
  // Aluminum - crystal end caps. This is where phonon hits are registered
  //
  G4VSolid* fAluminumSolid = new G4Tubs("aluminiumSolid",0.*cm,3.81*cm,0.01*cm,
                                        0.*deg, 360.*deg);

  G4LogicalVolume* fAluminumLogical =
    new G4LogicalVolume(fAluminumSolid,fAluminum,"fAluminumLogical");
  G4VPhysicalVolume* aluminumTopPhysical = new G4PVPlacement(0,
    G4ThreeVector(0.,0.,1.28*cm), fAluminumLogical, "fAluminumPhysical",
    worldLogical,false,0);
  G4VPhysicalVolume* aluminumBotPhysical = new G4PVPlacement(0,
    G4ThreeVector(0.,0.,-1.28*cm), fAluminumLogical, "fAluminumPhysical",
    worldLogical,false,1);

  //
  // detector -- Note : Aluminum electrode sensitivity is attached to Germanium 
  //
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  if (!electrodeSensitivity)
    electrodeSensitivity = new PhononSensitivity("PhononElectrode");
  SDman->AddNewDetector(electrodeSensitivity);
  fGermaniumLogical->SetSensitiveDetector(electrodeSensitivity);

  //
  // surface between Al and Ge determines phonon reflection/absorption

  if (!fConstructed) {
    const G4double GHz = 1e9 * hertz; 

    //the following coefficients and cutoff values are not well-motivated
    //the code below is used only to demonstrate how to set these values.
    const std::vector<G4double> anhCoeffs = {0, 0, 0, 0, 0, 1.51e-14};
    const std::vector<G4double> diffCoeffs =
      {5.88e-2, 7.83e-4, -2.47e-6, 1.71e-8, -2.98e-11};
    const std::vector<G4double> specCoeffs =
      {0,928, -2.03e-4, -3.21e-6, 3.1e-9, 2.9e-13};

    const G4double anhCutoff = 520., reflCutoff = 350.;   // Units external

    topSurfProp = new G4CMPSurfaceProperty("TopAlSurf", 1.0, 0.0, 0.0, 0.0,
					  	        0.3, 1.0, 0.0, 0.0);
    topSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					 diffCoeffs, specCoeffs, GHz, GHz, GHz);

    botSurfProp = new G4CMPSurfaceProperty("BotAlSurf", 1.0, 0.0, 0.0, 0.0,
					   	        0.3, 1.0, 0.0, 0.0);
    botSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					 diffCoeffs, specCoeffs, GHz, GHz, GHz);

    wallSurfProp = new G4CMPSurfaceProperty("WallAlSurf", 0.0, 1.0, 0.0, 0.0,
					    	          0.0, 1.0, 0.0, 0.0);
    wallSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					  diffCoeffs, specCoeffs, GHz, GHz,GHz);
  }

  new G4CMPLogicalBorderSurface("iZIPTop", GePhys, aluminumTopPhysical,
				topSurfProp);
  new G4CMPLogicalBorderSurface("iZIPBot", GePhys, aluminumBotPhysical,
				botSurfProp);
  new G4CMPLogicalBorderSurface("iZIPWall", GePhys, fWorldPhys,
				wallSurfProp);

  //                                        
  // Visualization attributes
  //
  worldLogical->SetVisAttributes(G4VisAttributes::Invisible);
  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  fGermaniumLogical->SetVisAttributes(simpleBoxVisAtt);
  fAluminumLogical->SetVisAttributes(simpleBoxVisAtt);
}

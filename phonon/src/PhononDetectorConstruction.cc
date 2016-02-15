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
/// \file exoticphysics/phonon/src/PhononDetectorConstruction.cc \brief
/// Implementation of the PhononDetectorConstruction class
//
// $Id: a2016d29cc7d1e75482bfc623a533d20b60390da $
//
// 20140321  Drop passing placement transform to G4LatticePhysical

#include "PhononDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4PVPlacement.hh"
#include "G4UniformMagField.hh"
#include "G4TransportationManager.hh"
//#include "G4GeometryManager.hh"
#include "G4FieldManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SDManager.hh"

#include "G4CMPElectrodeSensitivity.hh"
#include "G4CMPSurfaceProperty.hh"
#include "G4LatticePhysical.hh"
#include "G4LatticeLogical.hh"
#include "G4LatticeManager.hh"

#include "G4UserLimits.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PhononDetectorConstruction::PhononDetectorConstruction()
 : fConstructed(false), fIfField(true) {
  fLiquidHelium = NULL;
  fGermanium = NULL;
  fAluminum = NULL;
  fTungsten = NULL;
  fWorldPhys = NULL;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PhononDetectorConstruction::~PhononDetectorConstruction() {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* PhononDetectorConstruction::Construct()
{
  if(!fConstructed)
  {
    //G4GeometryManager::GetInstance()->SetWorldMaximumExtent(1e-6*mm);
    fConstructed = true;
    DefineMaterials();
    SetupGeometry();
  }
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
  GePhysical->SetMillerOrientation(0,0,1);
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
  G4CMPElectrodeSensitivity* electrodeSensitivity =
    new G4CMPElectrodeSensitivity("G4CMPElectrode");
  SDman->AddNewDetector(electrodeSensitivity);
  fGermaniumLogical->SetSensitiveDetector(electrodeSensitivity);

  //
  // surface between Al and Ge determines phonon reflection/absorption
  //
  G4CMPSurfaceProperty* topSurfProp = new G4CMPSurfaceProperty("TopAlSurf", 0.30,
    0.0, 347.43e-6*eV, 3., 242e-12*s, 3.26e3*km/s, 0.29, 300e-9*m);
  G4CMPSurfaceProperty* botSurfProp = new G4CMPSurfaceProperty("BotAlSurf", 0.30,
    0.0, 347.43e-6*eV, 3., 242e-12*s, 3.26e3*km/s, 0.29, 300e-9*m);
  G4CMPSurfaceProperty* wallSurfProp = new G4CMPSurfaceProperty("WallAlSurf", 0.0,
    0.0, 347.43e-6*eV, 0.0, 0.0, 0.0, 0.0, 0.0);

  new G4LogicalBorderSurface("iZIPTop", GePhys, aluminumTopPhysical,
                                        topSurfProp);
  new G4LogicalBorderSurface("iZIPBot", GePhys, aluminumBotPhysical,
                                        botSurfProp);
  new G4LogicalBorderSurface("iZIPWall", GePhys, fWorldPhys,
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

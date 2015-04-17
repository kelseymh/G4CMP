// $Id$

#include "ChargeDetectorConstruction.hh"
#include "ChargeEMField.hh"
#include "G4Box.hh"
#include "G4IntersectionSolid.hh"
#include "ChargeElectrodeSensitivity.hh"
#include "G4Colour.hh"
#include "G4LatticeLogical.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"


ChargeDetectorConstruction::ChargeDetectorConstruction()
   : liquidHelium(0), germanium(0), aluminum(0), tungsten(0), worldPhys(0),
     constructed(false), ifField(true),
     latManager(G4LatticeManager::GetLatticeManager()) {;}

ChargeDetectorConstruction::~ChargeDetectorConstruction() {;}

G4VPhysicalVolume* ChargeDetectorConstruction::Construct() {
  if (!constructed) { 
    DefineMaterials();
    SetupGeometry();
    constructed = true;
  }
  return worldPhys;
}

void ChargeDetectorConstruction::DefineMaterials() { 
  G4NistManager* nistManager = G4NistManager::Instance();

  liquidHelium = nistManager->FindOrBuildMaterial("G4_AIR"); //FIXME
  germanium = nistManager->FindOrBuildMaterial("G4_Ge");
  aluminum = nistManager->FindOrBuildMaterial("G4_Al");
  tungsten = nistManager->FindOrBuildMaterial("G4_W");

  // Attach lattice information for germanium
  latManager->LoadLattice(germanium, "Ge");
}

void ChargeDetectorConstruction::SetupGeometry()
{
  // World
  G4VSolid* worldSolid = new G4Box("World",16.*cm,16.*cm,16.*cm);
  G4LogicalVolume* worldLogical = new G4LogicalVolume(worldSolid,
                                                      liquidHelium, "World");
  worldPhys = new G4PVPlacement(0, G4ThreeVector(), worldLogical, "World", 0,
                                false, 0);

  // Germanium
  G4VSolid* germaniumSolid = new G4Tubs("germaniumCyl", 0.*cm, 3.81*cm,
                                        1.27*cm, 0.*deg, 360.*deg);
  G4VSolid* zipCutBox = new G4Box("ZipCutBox", 37.7444*mm, 36.0934*mm,
                                  1.28*cm);
  G4VSolid* zipSolid = new G4IntersectionSolid("germaniumSolid", germaniumSolid,
                                               zipCutBox);
  G4LogicalVolume* germaniumLogical = new G4LogicalVolume(zipSolid, germanium,
                                                          "germaniumLogical");
  G4VPhysicalVolume* germaniumPhysical = new G4PVPlacement(0, G4ThreeVector(),
                                                           germaniumLogical,
                                                           "germaniumPhysical",
                                                           worldLogical,
                                                           false, 0);

  // Attach E field to germanium (logical volume, so all placements)
  // NOTE:  Pointer isn't used here; must be permanent to support macros
  if (ifField) new ChargeEMField(germaniumLogical);

  // Physical lattice for each placed detector
  G4LatticePhysical* detLattice =
    new G4LatticePhysical(latManager->GetLattice(germanium));
  detLattice->SetLatticeOrientation(0.,45.*deg);	// Flats at [110]
  latManager->RegisterLattice(germaniumPhysical, detLattice);

  // Aluminum
  G4VSolid* aluminumSolid = new G4Tubs("aluminiumSolid", 0.*cm, 3.81*cm,
                                       0.01*cm, 0.*deg, 360.*deg);
  G4LogicalVolume* aluminumLogical = new G4LogicalVolume(aluminumSolid,
                                                         aluminum,
                                                         "aluminumLogical");
  new G4PVPlacement(0, G4ThreeVector(0.,0.,1.28*cm), aluminumLogical,
                    "aluminumPhysical", worldLogical, false, 0);
  new G4PVPlacement(0, G4ThreeVector(0.,0.,-1.28*cm), aluminumLogical,
                    "aluminumPhysical", worldLogical, false, 1);

  // detector -- Note : Aluminum electrode sensitivity is attached to Germanium 
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  ChargeElectrodeSensitivity* electrodeSensitivity =
                              new ChargeElectrodeSensitivity("ChargeElectrode");
  SDman->AddNewDetector(electrodeSensitivity);
  germaniumLogical->SetSensitiveDetector(electrodeSensitivity);

  // Visualization attributes
  worldLogical->SetVisAttributes(G4VisAttributes::Invisible);
  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  germaniumLogical->SetVisAttributes(simpleBoxVisAtt);
  aluminumLogical->SetVisAttributes(simpleBoxVisAtt);
}



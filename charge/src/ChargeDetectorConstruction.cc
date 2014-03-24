// $Id$

#include "ChargeDetectorConstruction.hh"
#include "ChargeEMField.hh"
#include "G4Box.hh"
#include "G4CMPElectrodeSensitivity.hh"
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
   : liquidHelium(0), germanium(0), alminum(0), tungsten(0), worldPhys(0),
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

  liquidHelium = nistManager->FindOrBuildMaterial("G4_AIR"); // to be corrected.......
  germanium = nistManager->FindOrBuildMaterial("G4_Ge");
  alminum = nistManager->FindOrBuildMaterial("G4_Al");
  tungsten = nistManager->FindOrBuildMaterial("G4_W");

  // Attach lattice information for germanium
  latManager->LoadLattice(germanium, "Ge");
}

void ChargeDetectorConstruction::SetupGeometry()
{
  //     
  // World
  //
  G4VSolid* worldSolid = new G4Box("World",16.*cm,16.*cm,16.*cm);
  G4LogicalVolume* worldLogical = new G4LogicalVolume(worldSolid,liquidHelium,"World");
  //worldLogical->SetUserLimits(new G4UserLimits(0.01*mm, DBL_MAX, DBL_MAX, 0, 0));
  worldPhys = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"World",0,false,0);

  //                               
  // Germanium
  //  
  G4VSolid* germaniumSolid = new G4Tubs("germaniumSolid",0.*cm,3.81*cm,1.27*cm, 0.*deg, 360.*deg);
  G4LogicalVolume* germaniumLogical = new G4LogicalVolume(germaniumSolid,germanium,"germaniumLogical");
  G4VPhysicalVolume* germaniumPhysical = new G4PVPlacement(0,G4ThreeVector(),germaniumLogical,"germaniumPhysical",worldLogical,false,0);

  // Attach E field to germanium
  if (ifField) ChargeEMField setField(germaniumLogical);

  // Physical lattice for placed detector
  G4LatticePhysical* detLattice =
    new G4LatticePhysical(latManager->GetLattice(germanium));
  detLattice->SetMillerOrientation(0,0,1);
  latManager->RegisterLattice(germaniumPhysical, detLattice);

  //
  // Aluminum
  //
  G4VSolid* alminumSolid = new G4Tubs("aluminiumSolid",0.*cm,3.81*cm,0.01*cm, 0.*deg, 360.*deg);
  G4LogicalVolume* alminumLogical = new G4LogicalVolume(alminumSolid,alminum,"alminumLogical");
  new G4PVPlacement(0,G4ThreeVector(0.,0.,1.28*cm),alminumLogical,"alminumPhysical",worldLogical,false,0);
  new G4PVPlacement(0,G4ThreeVector(0.,0.,-1.28*cm),alminumLogical,"alminumPhysical",worldLogical,false,1);

  //
  // Tungsten
  //
  /*  G4VSolid* tungstenSolid = new G4Box("tungstenSolid",5.*mm,1.*mm,1.*mm);
  G4LogicalVolume* tungstenLogical = new G4LogicalVolume(tungstenSolid,tungsten,"tungstenLogical");
  new G4PVPlacement(0,G4ThreeVector(0.,0.6*cm,1.37*cm),tungstenLogical,"tungstenPhysical",worldLogical,false,0);
  new G4PVPlacement(0,G4ThreeVector(0.,0.6*cm,-1.37*cm),tungstenLogical,"tungstenPhysical",worldLogical,false,1);
  */

  //
  // detector -- Note : Aluminum electrode sensitivity is attached to Germanium 
  //
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4CMPElectrodeSensitivity* electrodeSensitivity = new G4CMPElectrodeSensitivity("G4CMPElectrode");
  SDman->AddNewDetector(electrodeSensitivity);
  germaniumLogical->SetSensitiveDetector(electrodeSensitivity);

  //                                        
  // Visualization attributes
  //
  worldLogical->SetVisAttributes(G4VisAttributes::Invisible);
  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  germaniumLogical->SetVisAttributes(simpleBoxVisAtt);
  alminumLogical->SetVisAttributes(simpleBoxVisAtt);
  //  tungstenLogical->SetVisAttributes(simpleBoxVisAtt);
}



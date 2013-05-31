
#include "Tst1DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UniformMagField.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SDManager.hh"

#include "G4ios.hh"

#include "AlminumElectrodeSensitivity.hh"
#include "PhysicalLattice.hh"
#include "LogicalLattice.hh"

#include "G4UserLimits.hh"

Tst1DetectorConstruction::Tst1DetectorConstruction():constructed(false),ifField(true)
{
  liquidHelium = NULL;
  germanium = NULL;
  alminum = NULL;
  tungsten = NULL;
  worldPhys = NULL;
	/*G4double position[4] = {0, 0, 0, -1.2*cm};
    G4double fieldVal[6];
    G4FieldManager* fieldMan = G4TransportationManager::GetTransportationManager()->GetFieldManager();
    const G4ElectricField* field = (G4ElectricField*)fieldMan->GetDetectorField();
    field->GetFieldValue(position,  fieldVal);
    G4cout <<  "Does field exist? : " <<  fieldMan->DoesFieldExist() <<  G4endl;
    G4cout <<  "Field Value:" << fieldVal[3] <<  " " <<  fieldVal[4] <<  " " <<  fieldVal[5] <<  G4endl;
    */
}

Tst1DetectorConstruction::~Tst1DetectorConstruction()
{
    delete field;
}



G4VPhysicalVolume* Tst1DetectorConstruction::Construct()
{
  if(!constructed)
  { 
    constructed = true;
    DefineMaterials();
    SetupGeometry();
  }
  return worldPhys;
}

void Tst1DetectorConstruction::DefineMaterials()
{ 
  G4NistManager* nistManager = G4NistManager::Instance();

  liquidHelium = nistManager->FindOrBuildMaterial("G4_AIR"); // to be corrected.......
  germanium = nistManager->FindOrBuildMaterial("G4_Ge");
  alminum = nistManager->FindOrBuildMaterial("G4_Al");
  tungsten = nistManager->FindOrBuildMaterial("G4_W");
}

void Tst1DetectorConstruction::SetupGeometry()
{
  //     
  // World
  //
  G4VSolid* worldSolid = new G4Box("World",16.*cm,16.*cm,16.*cm);
  G4LogicalVolume* worldLogical = new G4LogicalVolume(worldSolid,liquidHelium,"World");
  //worldLogical->SetUserLimits(new G4UserLimits(0.01*mm, DBL_MAX, DBL_MAX, 0, 0));
  worldPhys = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"World",0,false,0);

  // E field
  field = new Tst1EMField(worldLogical);

  
  //                               
  // Germanium
  //  
  G4VSolid* germaniumSolid = new G4Tubs("germaniumSolid",0.*cm,3.81*cm,1.27*cm, 0.*deg, 360.*deg);
  G4LogicalVolume* germaniumLogical = new G4LogicalVolume(germaniumSolid,germanium,"germaniumLogical");
  G4VPhysicalVolume* GePhys = new G4PVPlacement(0,G4ThreeVector(),germaniumLogical,"germaniumPhysical",worldLogical,false,0);

  //  G4cout<<"\nUserlimits:"<<germaniumLogical->GetUserLimits();

  //germaniumLogical->SetUserLimits(new G4UserLimits(1e-2*mm, DBL_MAX, DBL_MAX));

  //G4cout<<"\nUserlimits:"<<germaniumLogical->GetUserLimits();

  /////////////////////////
  //Germanium lattice information
  ////////////////////////
  /*
  LogicalLattice GeLogical;
 //Convention for polarization state: 0=LON, 1=ST, 2=FT
  if(GeLogical.load_NMap(161, 321, 0, "./CrystalMaps/LVec.ssv")){
    if(GeLogical.load_NMap(161, 321, 1, "./CrystalMaps/STVec.ssv")){
      if(GeLogical.load_NMap(161, 321, 2, "./CrystalMaps/FTVec.ssv")){
	G4cout<<"\nTst1DetectorConstruction::Loaded all three maps";}}}

  if(GeLogical.loadMap(161, 321, 0, "./CrystalMaps/L.ssv")){
    if(GeLogical.loadMap(161, 321, 1, "./CrystalMaps/ST.ssv")){
      if(GeLogical.loadMap(161, 321, 2, "./CrystalMaps/FT.ssv")){
	G4cout<<"\nTst1DetectorConstruction::Loaded all three velocity maps";}}}

  GeLogical.setDynamicalConstants(-0.732, -0.708, 0.376, 0.561);
  GeLogical.setScatteringConstant(3.67e-41*s*s*s);
  GeLogical.setAnhDecConstant(1.6456e-54*s*s*s*s);
  GeLogical.setLDOS(0.097834);
  GeLogical.setSTDOS(0.53539);
  GeLogical.setFTDOS(0.36677);
  
  PhysicalLattice GePhysical(GePhys,&GeLogical);
  //GePhysical.setMillerOrientation(0,0,1);
  //GePhysical.setLatticeOrientation((3.1/4)*rad,45);
  LatticeManager2::registerLattice(&GePhysical);
  */
  //
  // Alminum
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
  // electric field
  //
  //if(ifField)
  //{
  //}

  //
  // detector -- Note : Alminum electrode sensitivity is attached to Germanium 
  //
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  AlminumElectrodeSensitivity* electrodeSensitivity = new AlminumElectrodeSensitivity("AlminumElectrode");
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



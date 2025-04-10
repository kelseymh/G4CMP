/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file exoticphysics/phonon/src/QuasiparticleDetectorConstruction.cc \brief
/// Implementation of the QuasiparticleDetectorConstruction class
//
// $Id: a2016d29cc7d1e75482bfc623a533d20b60390da $
//
// 20140321  Drop passing placement transform to G4LatticePhysical
// 20211207  Replace G4Logical*Surface with G4CMP-specific versions.
// 20220809  [ For M. Hui ] -- Add frequency dependent surface properties.
// 20221006  Remove unused features; add phonon sensor pad with use of
//		G4CMPPhononElectrode to demonstrate KaplanQP.

#include "QuasiparticleDetectorConstruction.hh"
#include "QuasiparticleDetectorParameters.hh"
#include "QuasiparticleSensitivity.hh"
#include "QuasiparticleQubitHousing.hh"
#include "QuasiparticleTransmissionLine.hh"
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


using namespace QuasiparticleDetectorParameters;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

QuasiparticleDetectorConstruction::QuasiparticleDetectorConstruction()
  : fLiquidHelium(0), fSilicon(0), fAluminum(0), fTungsten(0), fCopper(0), fNiobium(0),
    fWorldPhys(0), topSurfProp(0), vacSurfProp(0), wallSurfProp(0), copTopSurfProp(0), copTopSurfProp2(0), topSurfProp2(0), alNbSurfProp(0),
    electrodeSensitivity(0), fConstructed(false) {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

QuasiparticleDetectorConstruction::~QuasiparticleDetectorConstruction() {
  delete topSurfProp;
  delete topSurfProp2;
  delete vacSurfProp;
  delete wallSurfProp;
  delete copTopSurfProp;
  delete copTopSurfProp2;
  delete alNbSurfProp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* QuasiparticleDetectorConstruction::Construct()
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

void QuasiparticleDetectorConstruction::DefineMaterials()
{ 
  G4NistManager* nistManager = G4NistManager::Instance();

  fLiquidHelium = nistManager->FindOrBuildMaterial("G4_AIR"); // to be corrected
  fSilicon = nistManager->FindOrBuildMaterial("G4_Si");
  fAluminum = nistManager->FindOrBuildMaterial("G4_Al");
  fCopper = nistManager->FindOrBuildMaterial("G4_Cu");
  fTungsten = nistManager->FindOrBuildMaterial("G4_W");
  fNiobium = nistManager->FindOrBuildMaterial("G4_Nb");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void QuasiparticleDetectorConstruction::SetupGeometry()
{
  //     
  // World
  //
  G4VSolid* worldSolid = new G4Box("World",16.*cm,16.*cm,16.*cm);
  G4LogicalVolume* worldLogical =
    new G4LogicalVolume(worldSolid,fLiquidHelium,"World");
  worldLogical->SetUserLimits(new G4UserLimits(100*mm, DBL_MAX, DBL_MAX, 0, 0));
  fWorldPhys = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"World",0,
                                 false,0);

  //A few different detector scenarios:
  //1. Geometry 1: Single Aluminum Cylinder in between Si and Cu
  //2. Geometry 2: Aluminum half-cylinder and Niobium half-cylinder in betwen Si and Cu
  //3. Geometry 3: Aluminum cylinder with concentric niobium sub-cylinder
  //4. Geometry 4: Geometry 3 with a rotated box offset from the center, within the Al
  //5. Geometry 5: A long thin strip of aluminum with a wider block of niobium
  //6. Geometry 6: Niobium cylinder with an offset aluminum cylinder, within which is a quarter-circle G4tubs of Niobium
  int geometryID = 7;
  
  if( geometryID == 1 ){
    
    //                               
    // Silicon cylinder - this is the volume in which we will propagate phonons
    //  
    G4VSolid* siliconSolid = new G4Tubs("siliconSolid",0.*cm,3.81*cm,1.27*cm, 0.*deg, 360.*deg);
    G4LogicalVolume* siliconLogical = new G4LogicalVolume(siliconSolid,fSilicon,"siliconLogical");
    G4VPhysicalVolume* SiPhys = new G4PVPlacement(0,G4ThreeVector(),siliconLogical,"siliconPhysical", worldLogical,false,0);
    
    
    //
    // Aluminum - crystal end caps. This is where phonon hits are registered
    //
    G4VSolid* aluminumSolid = new G4Tubs("aluminiumSolid",0.*cm,3.81*cm,4000000*nm,0.*deg, 360.*deg);
    G4LogicalVolume* aluminumLogical = new G4LogicalVolume(aluminumSolid,fAluminum,"aluminumLogical");
    G4VPhysicalVolume* aluminumTopPhysical = new G4PVPlacement(0,G4ThreeVector(0.,0.,(1.27*cm + 4000000 * nm)), aluminumLogical,
							       "aluminumPhysical",worldLogical,false,0);
  
    //
    // Copper, set to the full circle and above the Nb/Al layers
    //  
    G4VSolid* copperSolid = new G4Tubs("aluminiumSolid",0.*cm,3.81*cm,0.01*cm,0.*deg, 360.*deg);
    G4LogicalVolume* copperLogical = new G4LogicalVolume(copperSolid,fCopper,"copperLogical");
    G4VPhysicalVolume* copperOverTopPhysical = new G4PVPlacement(0,G4ThreeVector(0.,0.,(1.27*cm + 8000000 * nm + 0.01 * cm)), copperLogical,
								 "copperLogical",worldLogical,false,0);


  
    //
    //Silicon lattice information
    //

    // G4LatticeManager gives physics processes access to lattices by volume
    G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();
    G4LatticeLogical* SiLogical = LM->LoadLattice(fSilicon, "Si");
    G4LatticeLogical* AlLogical = LM->LoadLattice(fAluminum, "Al");
    G4LatticeLogical* CuLogical = LM->LoadLattice(fCopper, "Cu");

    // G4LatticePhysical assigns G4LatticeLogical a physical orientation
    G4LatticePhysical* SiPhysical = new G4LatticePhysical(SiLogical);
    G4LatticePhysical* AlPhysical = new G4LatticePhysical(AlLogical);
    G4LatticePhysical* CuPhysical = new G4LatticePhysical(CuLogical);  
    SiPhysical->SetMillerOrientation(1,0,0);
    AlPhysical->SetMillerOrientation(1,0,0);
    CuPhysical->SetMillerOrientation(1,0,0);
    LM->RegisterLattice(SiPhys, SiPhysical);
    LM->RegisterLattice(aluminumTopPhysical, AlPhysical);
    LM->RegisterLattice(copperOverTopPhysical, CuPhysical);

    G4cout << "Al lattice LDOS: " << AlPhysical->GetLDOS() << ", scat rate: " << AlPhysical->GetPolycrystalElasticScatterMFP() << G4endl;

    //
    // surface between Al and Si determines phonon reflection/absorption
    //
    if (!fConstructed) {
      const G4double GHz = 1e9 * hertz; 

      //the following coefficients and cutoff values are not well-motivated
      //the code below is used only to demonstrate how to set these values.
      const std::vector<G4double> anhCoeffs = {0, 0, 0, 0, 0, 1.51e-14};
      const std::vector<G4double> diffCoeffs = {};
      const std::vector<G4double> specCoeffs = {};
    
      const G4double anhCutoff = 520., reflCutoff = 350.;   // Units external
      vacSurfProp = new G4CMPSurfaceProperty("VacCopSurf",
					     0.0, 1.0, 0.0, 0.0,
					     0.0, 1.0, 0.0, 0.0,
					     0.0, 1.0);
      vacSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					   diffCoeffs, specCoeffs, GHz, GHz, GHz);

      copTopSurfProp = new G4CMPSurfaceProperty("CopAlSurf",
						0.0, 1.0, 0.0, 0.0,
						0.0, 0.0, 0.0, 0.0,
						0.0, 1.0);
      copTopSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					      diffCoeffs, specCoeffs, GHz, GHz, GHz);
    
    
      topSurfProp = new G4CMPSurfaceProperty("AlSiSurf",
					     0.0, 1.0, 0.0, 0.0,
					     0.0, 0.0, 0.0, 0.0,
					     0.0, 1.0);
      topSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					   diffCoeffs, specCoeffs, GHz, GHz, GHz);
    
      wallSurfProp = new G4CMPSurfaceProperty("WallSurf",
					      0.0, 1.0, 0.0, 0.0,
					      0.0, 1.0, 0.0, 0.0,
					      0.0, 1.0);
      wallSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					    diffCoeffs, specCoeffs, GHz, GHz,GHz);

    }


    //
    // Separate surfaces for sensors vs. bare sidewall
    //
    //new G4CMPLogicalBorderSurface("AlSi", SiPhys, aluminumTopPhysical, vacSurfProp);
    new G4CMPLogicalBorderSurface("AlSi", SiPhys, aluminumTopPhysical,topSurfProp);
    new G4CMPLogicalBorderSurface("SiAl", aluminumTopPhysical, SiPhys, topSurfProp);
    //new G4CMPLogicalBorderSurface("SiAl", aluminumTopPhysical, SiPhys, vacSurfProp);

    new G4CMPLogicalBorderSurface("CuAl", copperOverTopPhysical, aluminumTopPhysical,copTopSurfProp);
    new G4CMPLogicalBorderSurface("AlCu", aluminumTopPhysical, copperOverTopPhysical,copTopSurfProp);

    new G4CMPLogicalBorderSurface("AlVac",aluminumTopPhysical,fWorldPhys,vacSurfProp);
    new G4CMPLogicalBorderSurface("VacAl",fWorldPhys,aluminumTopPhysical,vacSurfProp);
    new G4CMPLogicalBorderSurface("CuVac",copperOverTopPhysical,fWorldPhys,vacSurfProp);
    new G4CMPLogicalBorderSurface("VacCu",fWorldPhys,copperOverTopPhysical,vacSurfProp);
  
    //new G4CMPLogicalBorderSurface("VacAl",fWorldPhys,aluminumTopPhysical,vacSurfProp); //With this volume we try to reflect twice at each interface. BAD?
  
    new G4CMPLogicalBorderSurface("VacSi",SiPhys,fWorldPhys,wallSurfProp);
    new G4CMPLogicalBorderSurface("SiVac",fWorldPhys,SiPhys,wallSurfProp);

  
    //                                        
    // Visualization attributes
    //
    worldLogical->SetVisAttributes(G4VisAttributes::Invisible);
    G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    simpleBoxVisAtt->SetVisibility(true);
    siliconLogical->SetVisAttributes(simpleBoxVisAtt);
    aluminumLogical->SetVisAttributes(simpleBoxVisAtt);
  }

  
  //-----------------------------------------------------------------------------------------------------------------------
  //Geometry 2:
  if( geometryID == 2 ){
    
    //                               
    // Silicon cylinder - this is the volume in which we will propagate phonons
    //  
    G4VSolid* siliconSolid = new G4Tubs("siliconSolid",0.*cm,3.81*cm,1.27*cm, 0.*deg, 360.*deg);
    G4LogicalVolume* siliconLogical = new G4LogicalVolume(siliconSolid,fSilicon,"siliconLogical");
    G4VPhysicalVolume* SiPhys = new G4PVPlacement(0,G4ThreeVector(),siliconLogical,"siliconPhysical", worldLogical,false,0);

  
    //
    // Aluminum - crystal end caps. This is where phonon hits are registered
    //
    G4VSolid* aluminumSolid = new G4Tubs("aluminiumSolid",0.*cm,3.81*cm,4000000*nm,0.*deg, 180.*deg);
    G4LogicalVolume* aluminumLogical = new G4LogicalVolume(aluminumSolid,fAluminum,"aluminumLogical");
    G4VPhysicalVolume* aluminumTopPhysical = new G4PVPlacement(0,G4ThreeVector(0.,0.,(1.27*cm + 4000000 * nm)), aluminumLogical,
							       "aluminumPhysical",worldLogical,false,0);
  
  
    //
    // Niobium - half circle at same z as aluminum (but with different delta)
    //
    G4VSolid* niobiumSolid = new G4Tubs("niobiumSolid",0.*cm,3.81*cm,4000000*nm,180.*deg, 180.*deg);  
    G4LogicalVolume* niobiumLogical = new G4LogicalVolume(niobiumSolid,fNiobium,"niobiumLogical");
    G4VPhysicalVolume* niobiumTopPhysical = new G4PVPlacement(0,G4ThreeVector(0.,0.,(1.27*cm + 4000000 * nm)),niobiumLogical,
							      "niobiumPhysical",worldLogical,false,0);
  
  
    //
    // Copper, set to the full circle and above the Nb/Al layers
    //  
    G4VSolid* copperSolid = new G4Tubs("aluminiumSolid",0.*cm,3.81*cm,0.01*cm,0.*deg, 360.*deg);
    G4LogicalVolume* copperLogical = new G4LogicalVolume(copperSolid,fCopper,"copperLogical");
    G4VPhysicalVolume* copperOverTopPhysical = new G4PVPlacement(0,G4ThreeVector(0.,0.,(1.27*cm + 8000000 * nm + 0.01 * cm)), copperLogical,
								 "copperLogical",worldLogical,false,0);


  
    //
    //Silicon lattice information
    //

    // G4LatticeManager gives physics processes access to lattices by volume
    G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();
    G4LatticeLogical* SiLogical = LM->LoadLattice(fSilicon, "Si");
    G4LatticeLogical* AlLogical = LM->LoadLattice(fAluminum, "Al");
    G4LatticeLogical* CuLogical = LM->LoadLattice(fCopper, "Cu");
    G4LatticeLogical* NbLogical = LM->LoadLattice(fNiobium, "Nb");
    
    // G4LatticePhysical assigns G4LatticeLogical a physical orientation
    G4LatticePhysical* SiPhysical = new G4LatticePhysical(SiLogical);
    G4LatticePhysical* AlPhysical = new G4LatticePhysical(AlLogical);
    G4LatticePhysical* CuPhysical = new G4LatticePhysical(CuLogical);
    G4LatticePhysical* NbPhysical = new G4LatticePhysical(NbLogical);
    SiPhysical->SetMillerOrientation(1,0,0);
    AlPhysical->SetMillerOrientation(1,0,0);
    CuPhysical->SetMillerOrientation(1,0,0);
    NbPhysical->SetMillerOrientation(1,0,0);
    LM->RegisterLattice(SiPhys, SiPhysical);
    LM->RegisterLattice(aluminumTopPhysical, AlPhysical);
    LM->RegisterLattice(copperOverTopPhysical, CuPhysical);
    LM->RegisterLattice(niobiumTopPhysical,NbPhysical);

    G4cout << "Al lattice LDOS: " << AlPhysical->GetLDOS() << ", scat rate: " << AlPhysical->GetPolycrystalElasticScatterMFP() << G4endl;

    //
    // surface between Al and Si determines phonon reflection/absorption
    //
    if (!fConstructed) {
      const G4double GHz = 1e9 * hertz; 

      //the following coefficients and cutoff values are not well-motivated
      //the code below is used only to demonstrate how to set these values.
      const std::vector<G4double> anhCoeffs = {0, 0, 0, 0, 0, 1.51e-14};
      const std::vector<G4double> diffCoeffs = {};
      const std::vector<G4double> specCoeffs = {};
      
    
      const G4double anhCutoff = 520., reflCutoff = 350.;   // Units external
      
      //For the the interface of the Cu and the world
      vacSurfProp = new G4CMPSurfaceProperty("VacCopSurf",
					     0.0, 1.0, 0.0, 0.0,
					     0.0, 1.0, 0.0, 0.0,
					     0.0, 1.0);
      vacSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					   diffCoeffs, specCoeffs, GHz, GHz, GHz);

      //For the the interface of the Al and Cu    
      copTopSurfProp = new G4CMPSurfaceProperty("CopAlSurf",
						0.0, 1.0, 0.0, 0.0,
						0.0, 0.0, 0.0, 0.0,
						0.0, 1.0);
      copTopSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					      diffCoeffs, specCoeffs, GHz, GHz, GHz);

      //For the the interface of the Nb and Cu    
      copTopSurfProp2 = new G4CMPSurfaceProperty("CopNbSurf",
						 0.0, 1.0, 0.0, 0.0,
						 0.0, 0.0, 0.0, 0.0,
						 0.0, 1.0);
      copTopSurfProp2->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					       diffCoeffs, specCoeffs, GHz, GHz, GHz);

      
      //For the the interface of the Al and Si    
      topSurfProp = new G4CMPSurfaceProperty("AlSiSurf",
					     0.0, 1.0, 0.0, 0.0,
					     1.0, 0.0, 0.0, 0.0,
					     0.0, 1.0);
      topSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					   diffCoeffs, specCoeffs, GHz, GHz, GHz);

      //For the the interface of the Nb and Si
      topSurfProp2 = new G4CMPSurfaceProperty("NbSiSurf",
					      0.0, 1.0, 0.0, 0.0,
					      1.0, 0.0, 0.0, 0.0,
					      0.0, 1.0);
      topSurfProp2->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					    diffCoeffs, specCoeffs, GHz, GHz, GHz);
      
      //For the aluminum/niobium interface
      alNbSurfProp = new G4CMPSurfaceProperty("AlNbSurf",
					      0.0,1.0,0.0,0.0,
					      0.0,0.0,0.0,0.0,
					      0.0,0.0);
      alNbSurfProp->AddScatteringProperties(anhCutoff,reflCutoff,anhCoeffs,
					    diffCoeffs,specCoeffs,GHz,GHz,GHz);
      
      //For the interface of anything with the world
      wallSurfProp = new G4CMPSurfaceProperty("WallSurf",
					      0.0, 1.0, 0.0, 0.0,
					      0.0, 1.0, 0.0, 0.0,
					      0.0, 1.0);
      wallSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					    diffCoeffs, specCoeffs, GHz, GHz,GHz);

      
      
    }


    //
    // Separate surfaces for sensors vs. bare sidewall
    //
    new G4CMPLogicalBorderSurface("AlSi", SiPhys, aluminumTopPhysical,topSurfProp);
    new G4CMPLogicalBorderSurface("SiAl", aluminumTopPhysical, SiPhys, topSurfProp);
    new G4CMPLogicalBorderSurface("NbSi", SiPhys, niobiumTopPhysical,topSurfProp2);
    new G4CMPLogicalBorderSurface("SiNb", niobiumTopPhysical, SiPhys, topSurfProp2);

    new G4CMPLogicalBorderSurface("NbAl", niobiumTopPhysical,aluminumTopPhysical,alNbSurfProp);
    new G4CMPLogicalBorderSurface("AlNb", aluminumTopPhysical,niobiumTopPhysical,alNbSurfProp);
    
    new G4CMPLogicalBorderSurface("CuAl", copperOverTopPhysical, aluminumTopPhysical,copTopSurfProp);
    new G4CMPLogicalBorderSurface("AlCu", aluminumTopPhysical, copperOverTopPhysical,copTopSurfProp);
    new G4CMPLogicalBorderSurface("CuNb", copperOverTopPhysical, niobiumTopPhysical,copTopSurfProp2);
    new G4CMPLogicalBorderSurface("NbCu", niobiumTopPhysical, copperOverTopPhysical,copTopSurfProp2);
    
    new G4CMPLogicalBorderSurface("AlVac",aluminumTopPhysical,fWorldPhys,vacSurfProp);
    new G4CMPLogicalBorderSurface("VacAl",fWorldPhys,aluminumTopPhysical,vacSurfProp);
    new G4CMPLogicalBorderSurface("NbVac",niobiumTopPhysical,fWorldPhys,vacSurfProp);
    new G4CMPLogicalBorderSurface("VacNb",fWorldPhys,niobiumTopPhysical,vacSurfProp);

    
    new G4CMPLogicalBorderSurface("CuVac",copperOverTopPhysical,fWorldPhys,vacSurfProp);
    new G4CMPLogicalBorderSurface("VacCu",fWorldPhys,copperOverTopPhysical,vacSurfProp);
  
    //new G4CMPLogicalBorderSurface("VacAl",fWorldPhys,aluminumTopPhysical,vacSurfProp); //With this volume we try to reflect twice at each interface. BAD?
  
    new G4CMPLogicalBorderSurface("VacSi",SiPhys,fWorldPhys,wallSurfProp);
    new G4CMPLogicalBorderSurface("SiVac",fWorldPhys,SiPhys,wallSurfProp);

  
    //                                        
    // Visualization attributes
    //
    worldLogical->SetVisAttributes(G4VisAttributes::Invisible);
    G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    simpleBoxVisAtt->SetVisibility(true);
    siliconLogical->SetVisAttributes(simpleBoxVisAtt);
    aluminumLogical->SetVisAttributes(simpleBoxVisAtt);
    niobiumLogical->SetVisAttributes(simpleBoxVisAtt);
  }

  //-----------------------------------------------------------------------------------------------------------------------
  //Geometry 3:
  if( geometryID == 3 ){
    
    //                               
    // Silicon cylinder - this is the volume in which we will propagate phonons
    //  
    G4VSolid* siliconSolid = new G4Tubs("siliconSolid",0.*cm,3.81*cm,1.27*cm, 0.*deg, 360.*deg);
    G4LogicalVolume* siliconLogical = new G4LogicalVolume(siliconSolid,fSilicon,"siliconLogical");
    G4VPhysicalVolume* SiPhys = new G4PVPlacement(0,G4ThreeVector(),siliconLogical,"siliconPhysical", worldLogical,false,0);

  
    //
    // Aluminum - crystal end caps. This is where phonon hits are registered
    //
    G4VSolid* aluminumSolid = new G4Tubs("aluminiumSolid",0.*cm,3.81*cm,4000000*nm,0.*deg, 360.*deg);
    G4LogicalVolume* aluminumLogical = new G4LogicalVolume(aluminumSolid,fAluminum,"aluminumLogical");
    G4VPhysicalVolume* aluminumTopPhysical = new G4PVPlacement(0,G4ThreeVector(0.,0.,(1.27*cm + 4000000 * nm)), aluminumLogical,
							       "aluminumPhysical",worldLogical,false,0);
  
  
    //
    // Niobium - full circle at same z as aluminum (but with different radius and delta)
    //
    G4VSolid* niobiumSolid = new G4Tubs("niobiumSolid",0.*cm,1.0*cm,4000000*nm,0.*deg, 360.*deg);  
    G4LogicalVolume* niobiumLogical = new G4LogicalVolume(niobiumSolid,fNiobium,"niobiumLogical");
    G4VPhysicalVolume* niobiumTopPhysical = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),niobiumLogical,
							      "niobiumPhysical",aluminumLogical,false,0);
  
  
    //
    // Copper, set to the full circle and above the Nb/Al layers
    //  
    G4VSolid* copperSolid = new G4Tubs("aluminiumSolid",0.*cm,3.81*cm,0.01*cm,0.*deg, 360.*deg);
    G4LogicalVolume* copperLogical = new G4LogicalVolume(copperSolid,fCopper,"copperLogical");
    G4VPhysicalVolume* copperOverTopPhysical = new G4PVPlacement(0,G4ThreeVector(0.,0.,(1.27*cm + 8000000 * nm + 0.01 * cm)), copperLogical,
								 "copperLogical",worldLogical,false,0);


  
    //
    //Silicon lattice information
    //

    // G4LatticeManager gives physics processes access to lattices by volume
    G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();
    G4LatticeLogical* SiLogical = LM->LoadLattice(fSilicon, "Si");
    G4LatticeLogical* AlLogical = LM->LoadLattice(fAluminum, "Al");
    G4LatticeLogical* CuLogical = LM->LoadLattice(fCopper, "Cu");
    G4LatticeLogical* NbLogical = LM->LoadLattice(fNiobium, "Nb");
    
    // G4LatticePhysical assigns G4LatticeLogical a physical orientation
    G4LatticePhysical* SiPhysical = new G4LatticePhysical(SiLogical);
    G4LatticePhysical* AlPhysical = new G4LatticePhysical(AlLogical);
    G4LatticePhysical* CuPhysical = new G4LatticePhysical(CuLogical);
    G4LatticePhysical* NbPhysical = new G4LatticePhysical(NbLogical);
    SiPhysical->SetMillerOrientation(1,0,0);
    AlPhysical->SetMillerOrientation(1,0,0);
    CuPhysical->SetMillerOrientation(1,0,0);
    NbPhysical->SetMillerOrientation(1,0,0);
    LM->RegisterLattice(SiPhys, SiPhysical);
    LM->RegisterLattice(aluminumTopPhysical, AlPhysical);
    LM->RegisterLattice(copperOverTopPhysical, CuPhysical);
    LM->RegisterLattice(niobiumTopPhysical,NbPhysical);

    G4cout << "Al lattice LDOS: " << AlPhysical->GetLDOS() << ", scat rate: " << AlPhysical->GetPolycrystalElasticScatterMFP() << G4endl;

    //
    // surface between Al and Si determines phonon reflection/absorption
    //
    if (!fConstructed) {
      const G4double GHz = 1e9 * hertz; 

      //the following coefficients and cutoff values are not well-motivated
      //the code below is used only to demonstrate how to set these values.
      const std::vector<G4double> anhCoeffs = {0, 0, 0, 0, 0, 1.51e-14};
      const std::vector<G4double> diffCoeffs = {};
      const std::vector<G4double> specCoeffs = {};
      
    
      const G4double anhCutoff = 520., reflCutoff = 350.;   // Units external
      
      //For the the interface of the Cu and the world
      vacSurfProp = new G4CMPSurfaceProperty("VacCopSurf",
					     0.0, 1.0, 0.0, 0.0,
					     0.0, 1.0, 0.0, 0.0,
					     0.0, 1.0);
      vacSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					   diffCoeffs, specCoeffs, GHz, GHz, GHz);

      //For the the interface of the Al and Cu    
      copTopSurfProp = new G4CMPSurfaceProperty("CopAlSurf",
						0.0, 1.0, 0.0, 0.0,
						0.0, 0.0, 0.0, 0.0,
						0.0, 1.0);
      copTopSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					      diffCoeffs, specCoeffs, GHz, GHz, GHz);

      //For the the interface of the Nb and Cu    
      copTopSurfProp2 = new G4CMPSurfaceProperty("CopNbSurf",
						 0.0, 1.0, 0.0, 0.0,
						 0.0, 0.0, 0.0, 0.0,
						 0.0, 1.0);
      copTopSurfProp2->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					       diffCoeffs, specCoeffs, GHz, GHz, GHz);

      
      //For the the interface of the Al and Si    
      topSurfProp = new G4CMPSurfaceProperty("AlSiSurf",
					     0.0, 1.0, 0.0, 0.0,
					     1.0, 0.0, 0.0, 0.0,
					     0.0, 1.0);
      topSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					   diffCoeffs, specCoeffs, GHz, GHz, GHz);

      //For the the interface of the Nb and Si
      topSurfProp2 = new G4CMPSurfaceProperty("NbSiSurf",
					      0.0, 1.0, 0.0, 0.0,
					      1.0, 0.0, 0.0, 0.0,
					      0.0, 1.0);
      topSurfProp2->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					    diffCoeffs, specCoeffs, GHz, GHz, GHz);
      
      //For the aluminum/niobium interface
      alNbSurfProp = new G4CMPSurfaceProperty("AlNbSurf",
					      0.0,1.0,0.0,0.0,
					      0.0,0.0,0.0,0.0,
					      0.0,0.0);
      alNbSurfProp->AddScatteringProperties(anhCutoff,reflCutoff,anhCoeffs,
					    diffCoeffs,specCoeffs,GHz,GHz,GHz);
      
      //For the interface of anything with the world
      wallSurfProp = new G4CMPSurfaceProperty("WallSurf",
					      0.0, 1.0, 0.0, 0.0,
					      0.0, 1.0, 0.0, 0.0,
					      0.0, 1.0);
      wallSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					    diffCoeffs, specCoeffs, GHz, GHz,GHz);

      
      
    }


    //
    // Separate surfaces for sensors vs. bare sidewall
    //
    new G4CMPLogicalBorderSurface("AlSi", SiPhys, aluminumTopPhysical,topSurfProp);
    new G4CMPLogicalBorderSurface("SiAl", aluminumTopPhysical, SiPhys, topSurfProp);
    new G4CMPLogicalBorderSurface("NbSi", SiPhys, niobiumTopPhysical,topSurfProp2);
    new G4CMPLogicalBorderSurface("SiNb", niobiumTopPhysical, SiPhys, topSurfProp2);

    new G4CMPLogicalBorderSurface("NbAl", niobiumTopPhysical,aluminumTopPhysical,alNbSurfProp);
    new G4CMPLogicalBorderSurface("AlNb", aluminumTopPhysical,niobiumTopPhysical,alNbSurfProp);
    
    new G4CMPLogicalBorderSurface("CuAl", copperOverTopPhysical, aluminumTopPhysical,copTopSurfProp);
    new G4CMPLogicalBorderSurface("AlCu", aluminumTopPhysical, copperOverTopPhysical,copTopSurfProp);
    new G4CMPLogicalBorderSurface("CuNb", copperOverTopPhysical, niobiumTopPhysical,copTopSurfProp2);
    new G4CMPLogicalBorderSurface("NbCu", niobiumTopPhysical, copperOverTopPhysical,copTopSurfProp2);
    
    new G4CMPLogicalBorderSurface("AlVac",aluminumTopPhysical,fWorldPhys,vacSurfProp);
    new G4CMPLogicalBorderSurface("VacAl",fWorldPhys,aluminumTopPhysical,vacSurfProp);
    //    new G4CMPLogicalBorderSurface("NbVac",niobiumTopPhysical,fWorldPhys,vacSurfProp);
    //new G4CMPLogicalBorderSurface("VacNb",fWorldPhys,niobiumTopPhysical,vacSurfProp);

    
    new G4CMPLogicalBorderSurface("CuVac",copperOverTopPhysical,fWorldPhys,vacSurfProp);
    new G4CMPLogicalBorderSurface("VacCu",fWorldPhys,copperOverTopPhysical,vacSurfProp);
  
    //new G4CMPLogicalBorderSurface("VacAl",fWorldPhys,aluminumTopPhysical,vacSurfProp); //With this volume we try to reflect twice at each interface. BAD?
  
    new G4CMPLogicalBorderSurface("VacSi",SiPhys,fWorldPhys,wallSurfProp);
    new G4CMPLogicalBorderSurface("SiVac",fWorldPhys,SiPhys,wallSurfProp);

  
    //                                        
    // Visualization attributes
    //
    worldLogical->SetVisAttributes(G4VisAttributes::Invisible);
    G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    simpleBoxVisAtt->SetVisibility(true);
    siliconLogical->SetVisAttributes(simpleBoxVisAtt);
    aluminumLogical->SetVisAttributes(simpleBoxVisAtt);
    niobiumLogical->SetVisAttributes(simpleBoxVisAtt);
  }



  //-----------------------------------------------------------------------------------------------------------------------
  //Geometry 4:
  if( geometryID == 4 ){
    
    //                               
    // Silicon cylinder - this is the volume in which we will propagate phonons
    //  
    G4VSolid* siliconSolid = new G4Tubs("siliconSolid",0.*cm,3.81*cm,1.27*cm, 0.*deg, 360.*deg);
    G4LogicalVolume* siliconLogical = new G4LogicalVolume(siliconSolid,fSilicon,"siliconLogical");
    G4VPhysicalVolume* SiPhys = new G4PVPlacement(0,G4ThreeVector(),siliconLogical,"siliconPhysical", worldLogical,false,0);

  
    //
    // Aluminum - crystal end caps. This is where phonon hits are registered
    //
    G4VSolid* aluminumSolid = new G4Tubs("aluminiumSolid",0.*cm,3.81*cm,4000000*nm,0.*deg, 360.*deg);
    G4LogicalVolume* aluminumLogical = new G4LogicalVolume(aluminumSolid,fAluminum,"aluminumLogical");
    G4VPhysicalVolume* aluminumTopPhysical = new G4PVPlacement(0,G4ThreeVector(0.,0.,(1.27*cm + 4000000 * nm)), aluminumLogical,
							       "aluminumPhysical",worldLogical,false,0);
  
  
    //
    // Niobium - full circle at same z as aluminum (but with different radius and delta)
    //
    G4VSolid* niobiumSolid = new G4Tubs("niobiumSolid",0.*cm,1.0*cm,4000000*nm,0.*deg, 360.*deg);  
    G4LogicalVolume* niobiumLogical = new G4LogicalVolume(niobiumSolid,fNiobium,"niobiumLogical");
    G4VPhysicalVolume* niobiumTopPhysical = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),niobiumLogical,
							      "niobiumPhysical",aluminumLogical,false,0);

    //
    // Niobium misc - box at same z as aluminum (offset from center and rotated around its own z axis)
    //
    G4RotationMatrix * rotLine = new G4RotationMatrix();
    rotLine->rotateZ(45*deg);
    G4VSolid* niobiumSolid2 = new G4Box("niobiumSolid2",0.5*cm,0.5*cm,4000000*nm);
    G4LogicalVolume* niobiumLogical2 = new G4LogicalVolume(niobiumSolid2,fNiobium,"niobiumLogical2");
    G4VPhysicalVolume* niobiumTopPhysical2 = new G4PVPlacement(rotLine,G4ThreeVector(1.9*cm,1.9*cm,0.),niobiumLogical2,
							      "niobiumPhysical2",aluminumLogical,false,0);
     
    //
    // Copper, set to the full circle and above the Nb/Al layers
    //
    G4VSolid* copperSolid = new G4Tubs("aluminiumSolid",0.*cm,3.81*cm,0.01*cm,0.*deg, 360.*deg);
    G4LogicalVolume* copperLogical = new G4LogicalVolume(copperSolid,fCopper,"copperLogical");
    G4VPhysicalVolume* copperOverTopPhysical = new G4PVPlacement(0,G4ThreeVector(0.,0.,(1.27*cm + 8000000 * nm + 0.01 * cm)), copperLogical,
								 "copperLogical",worldLogical,false,0);


  
    //
    //Silicon lattice information
    //

    // G4LatticeManager gives physics processes access to lattices by volume
    G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();
    G4LatticeLogical* SiLogical = LM->LoadLattice(fSilicon, "Si");
    G4LatticeLogical* AlLogical = LM->LoadLattice(fAluminum, "Al");
    G4LatticeLogical* CuLogical = LM->LoadLattice(fCopper, "Cu");
    G4LatticeLogical* NbLogical = LM->LoadLattice(fNiobium, "Nb");
    G4LatticeLogical* NbLogical2 = LM->LoadLattice(fNiobium, "Nb");
    
    // G4LatticePhysical assigns G4LatticeLogical a physical orientation
    G4LatticePhysical* SiPhysical = new G4LatticePhysical(SiLogical);
    G4LatticePhysical* AlPhysical = new G4LatticePhysical(AlLogical);
    G4LatticePhysical* CuPhysical = new G4LatticePhysical(CuLogical);
    G4LatticePhysical* NbPhysical = new G4LatticePhysical(NbLogical);
    G4LatticePhysical* NbPhysical2 = new G4LatticePhysical(NbLogical2);
    SiPhysical->SetMillerOrientation(1,0,0);
    AlPhysical->SetMillerOrientation(1,0,0);
    CuPhysical->SetMillerOrientation(1,0,0);
    NbPhysical->SetMillerOrientation(1,0,0);
    NbPhysical2->SetMillerOrientation(1,0,0);
    LM->RegisterLattice(SiPhys, SiPhysical);
    LM->RegisterLattice(aluminumTopPhysical, AlPhysical);
    LM->RegisterLattice(copperOverTopPhysical, CuPhysical);
    LM->RegisterLattice(niobiumTopPhysical,NbPhysical);
    LM->RegisterLattice(niobiumTopPhysical2,NbPhysical2);

    G4cout << "Al lattice LDOS: " << AlPhysical->GetLDOS() << ", scat rate: " << AlPhysical->GetPolycrystalElasticScatterMFP() << G4endl;

    //
    // surface between Al and Si determines phonon reflection/absorption
    //
    if (!fConstructed) {
      const G4double GHz = 1e9 * hertz; 

      //the following coefficients and cutoff values are not well-motivated
      //the code below is used only to demonstrate how to set these values.
      const std::vector<G4double> anhCoeffs = {0, 0, 0, 0, 0, 1.51e-14};
      const std::vector<G4double> diffCoeffs = {};
      const std::vector<G4double> specCoeffs = {};
      
    
      const G4double anhCutoff = 520., reflCutoff = 350.;   // Units external
      
      //For the the interface of the Cu and the world
      vacSurfProp = new G4CMPSurfaceProperty("VacCopSurf",
					     0.0, 1.0, 0.0, 0.0,
					     0.0, 1.0, 0.0, 0.0,
					     0.0, 1.0);
      vacSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					   diffCoeffs, specCoeffs, GHz, GHz, GHz);

      //For the the interface of the Al and Cu    
      copTopSurfProp = new G4CMPSurfaceProperty("CopAlSurf",
						0.0, 1.0, 0.0, 0.0,
						0.0, 0.0, 0.0, 0.0,
						0.0, 1.0);
      copTopSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					      diffCoeffs, specCoeffs, GHz, GHz, GHz);

      //For the the interface of the Nb and Cu    
      copTopSurfProp2 = new G4CMPSurfaceProperty("CopNbSurf",
						 0.0, 1.0, 0.0, 0.0,
						 0.0, 0.0, 0.0, 0.0,
						 0.0, 1.0);
      copTopSurfProp2->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					       diffCoeffs, specCoeffs, GHz, GHz, GHz);

      
      //For the the interface of the Al and Si    
      topSurfProp = new G4CMPSurfaceProperty("AlSiSurf",
					     0.0, 1.0, 0.0, 0.0,
					     1.0, 0.0, 0.0, 0.0,
					     0.0, 1.0);
      topSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					   diffCoeffs, specCoeffs, GHz, GHz, GHz);

      //For the the interface of the Nb and Si
      topSurfProp2 = new G4CMPSurfaceProperty("NbSiSurf",
					      0.0, 1.0, 0.0, 0.0,
					      1.0, 0.0, 0.0, 0.0,
					      0.0, 1.0);
      topSurfProp2->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					    diffCoeffs, specCoeffs, GHz, GHz, GHz);
      
      //For the aluminum/niobium interface
      alNbSurfProp = new G4CMPSurfaceProperty("AlNbSurf",
					      0.0,1.0,0.0,0.0,
					      0.0,0.0,0.0,0.0,
					      0.0,0.0);
      alNbSurfProp->AddScatteringProperties(anhCutoff,reflCutoff,anhCoeffs,
					    diffCoeffs,specCoeffs,GHz,GHz,GHz);
      
      //For the interface of anything with the world
      wallSurfProp = new G4CMPSurfaceProperty("WallSurf",
					      0.0, 1.0, 0.0, 0.0,
					      0.0, 1.0, 0.0, 0.0,
					      0.0, 1.0);
      wallSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					    diffCoeffs, specCoeffs, GHz, GHz,GHz);

      
      
    }


    //
    // Separate surfaces for sensors vs. bare sidewall
    //
    new G4CMPLogicalBorderSurface("AlSi", SiPhys, aluminumTopPhysical,topSurfProp);
    new G4CMPLogicalBorderSurface("SiAl", aluminumTopPhysical, SiPhys, topSurfProp);
    new G4CMPLogicalBorderSurface("NbSi", SiPhys, niobiumTopPhysical,topSurfProp2);
    new G4CMPLogicalBorderSurface("NbSi2", SiPhys, niobiumTopPhysical2,topSurfProp2);
    new G4CMPLogicalBorderSurface("SiNb", niobiumTopPhysical, SiPhys, topSurfProp2);
    new G4CMPLogicalBorderSurface("SiNb2", niobiumTopPhysical2, SiPhys, topSurfProp2);

    new G4CMPLogicalBorderSurface("NbAl", niobiumTopPhysical,aluminumTopPhysical,alNbSurfProp);
    new G4CMPLogicalBorderSurface("AlNb", aluminumTopPhysical,niobiumTopPhysical,alNbSurfProp);
    new G4CMPLogicalBorderSurface("NbAl2", niobiumTopPhysical2,aluminumTopPhysical,alNbSurfProp);
    new G4CMPLogicalBorderSurface("AlNb2", aluminumTopPhysical,niobiumTopPhysical2,alNbSurfProp);


    
    new G4CMPLogicalBorderSurface("CuAl", copperOverTopPhysical, aluminumTopPhysical,copTopSurfProp);
    new G4CMPLogicalBorderSurface("AlCu", aluminumTopPhysical, copperOverTopPhysical,copTopSurfProp);
    new G4CMPLogicalBorderSurface("CuNb", copperOverTopPhysical, niobiumTopPhysical,copTopSurfProp2);
    new G4CMPLogicalBorderSurface("CuNb2", copperOverTopPhysical, niobiumTopPhysical2,copTopSurfProp2);
    new G4CMPLogicalBorderSurface("NbCu", niobiumTopPhysical, copperOverTopPhysical,copTopSurfProp2);
    new G4CMPLogicalBorderSurface("NbCu2", niobiumTopPhysical2, copperOverTopPhysical,copTopSurfProp2);
    
    new G4CMPLogicalBorderSurface("AlVac",aluminumTopPhysical,fWorldPhys,vacSurfProp);
    new G4CMPLogicalBorderSurface("VacAl",fWorldPhys,aluminumTopPhysical,vacSurfProp);
    
    new G4CMPLogicalBorderSurface("CuVac",copperOverTopPhysical,fWorldPhys,vacSurfProp);
    new G4CMPLogicalBorderSurface("VacCu",fWorldPhys,copperOverTopPhysical,vacSurfProp);
  
    //new G4CMPLogicalBorderSurface("VacAl",fWorldPhys,aluminumTopPhysical,vacSurfProp); //With this volume we try to reflect twice at each interface. BAD?
  
    new G4CMPLogicalBorderSurface("VacSi",SiPhys,fWorldPhys,wallSurfProp);
    new G4CMPLogicalBorderSurface("SiVac",fWorldPhys,SiPhys,wallSurfProp);

  
    //                                        
    // Visualization attributes
    //
    worldLogical->SetVisAttributes(G4VisAttributes::Invisible);
    G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    simpleBoxVisAtt->SetVisibility(true);
    siliconLogical->SetVisAttributes(simpleBoxVisAtt);
    aluminumLogical->SetVisAttributes(simpleBoxVisAtt);
    niobiumLogical->SetVisAttributes(simpleBoxVisAtt);
    niobiumLogical2->SetVisAttributes(simpleBoxVisAtt);
  }
  


  //-----------------------------------------------------------------------------------------------------------------------
  //Geometry 5:
  if( geometryID == 5 ){
    
    //                               
    // Silicon cylinder - this is the volume in which we will propagate phonons
    //  
    G4VSolid* siliconSolid = new G4Tubs("siliconSolid",0.*cm,3.81*cm,1.27*cm, 0.*deg, 360.*deg);
    G4LogicalVolume* siliconLogical = new G4LogicalVolume(siliconSolid,fSilicon,"siliconLogical");
    G4VPhysicalVolume* SiPhys = new G4PVPlacement(0,G4ThreeVector(),siliconLogical,"siliconPhysical", worldLogical,false,0);

  
    //
    // Aluminum - thin rectangular slab
    //
    G4VSolid* aluminumSolid = new G4Box("aluminumSolid",3.0*cm,0.002*cm,4000000*nm);
    G4LogicalVolume* aluminumLogical = new G4LogicalVolume(aluminumSolid,fAluminum,"aluminumLogical");
    G4VPhysicalVolume* aluminumTopPhysical = new G4PVPlacement(0,G4ThreeVector(0.,0.,(1.27*cm + 4000000 * nm)), aluminumLogical,
							       "aluminumPhysical",worldLogical,false,0);

    //
    // Niobium - thin rectangular slab on end of aluminum
    //
    G4VSolid* niobiumSolid = new G4Box("niobiumSolid",0.5*cm,0.2*cm,4000000*nm);
    G4LogicalVolume* niobiumLogical = new G4LogicalVolume(niobiumSolid,fNiobium,"niobiumLogical");
    G4VPhysicalVolume* niobiumTopPhysical = new G4PVPlacement(0,G4ThreeVector(3.5*cm,0.0,(1.27*cm + 4000000 * nm)), niobiumLogical,
							       "niobiumPhysical",worldLogical,false,0);
  
  
    //
    // Copper, set to the full circle and above the Nb/Al layers
    //
    G4VSolid* copperSolid = new G4Tubs("aluminiumSolid",0.*cm,3.81*cm,0.01*cm,0.*deg, 360.*deg);
    G4LogicalVolume* copperLogical = new G4LogicalVolume(copperSolid,fCopper,"copperLogical");
    G4VPhysicalVolume* copperOverTopPhysical = new G4PVPlacement(0,G4ThreeVector(0.,0.,(1.27*cm + 8000000 * nm + 0.01 * cm)), copperLogical,
								 "copperLogical",worldLogical,false,0);


  
    //
    //Silicon lattice information
    //

    // G4LatticeManager gives physics processes access to lattices by volume
    G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();
    G4LatticeLogical* SiLogical = LM->LoadLattice(fSilicon, "Si");
    G4LatticeLogical* AlLogical = LM->LoadLattice(fAluminum, "Al");
    G4LatticeLogical* CuLogical = LM->LoadLattice(fCopper, "Cu");
    G4LatticeLogical* NbLogical = LM->LoadLattice(fNiobium, "Nb");
    G4LatticeLogical* NbLogical2 = LM->LoadLattice(fNiobium, "Nb");
    
    // G4LatticePhysical assigns G4LatticeLogical a physical orientation
    G4LatticePhysical* SiPhysical = new G4LatticePhysical(SiLogical);
    G4LatticePhysical* AlPhysical = new G4LatticePhysical(AlLogical);
    G4LatticePhysical* CuPhysical = new G4LatticePhysical(CuLogical);
    G4LatticePhysical* NbPhysical = new G4LatticePhysical(NbLogical);
    G4LatticePhysical* NbPhysical2 = new G4LatticePhysical(NbLogical2);
    SiPhysical->SetMillerOrientation(1,0,0);
    AlPhysical->SetMillerOrientation(1,0,0);
    CuPhysical->SetMillerOrientation(1,0,0);
    NbPhysical->SetMillerOrientation(1,0,0);
    NbPhysical2->SetMillerOrientation(1,0,0);
    LM->RegisterLattice(SiPhys, SiPhysical);
    LM->RegisterLattice(aluminumTopPhysical, AlPhysical);
    LM->RegisterLattice(copperOverTopPhysical, CuPhysical);
    LM->RegisterLattice(niobiumTopPhysical,NbPhysical);

    G4cout << "Al lattice LDOS: " << AlPhysical->GetLDOS() << ", scat rate: " << AlPhysical->GetPolycrystalElasticScatterMFP() << G4endl;

    //
    // surface between Al and Si determines phonon reflection/absorption
    //
    if (!fConstructed) {
      const G4double GHz = 1e9 * hertz; 

      //the following coefficients and cutoff values are not well-motivated
      //the code below is used only to demonstrate how to set these values.
      const std::vector<G4double> anhCoeffs = {0, 0, 0, 0, 0, 1.51e-14};
      const std::vector<G4double> diffCoeffs = {};
      const std::vector<G4double> specCoeffs = {};
      
    
      const G4double anhCutoff = 520., reflCutoff = 350.;   // Units external
      
      //For the the interface of the Cu and the world
      vacSurfProp = new G4CMPSurfaceProperty("VacCopSurf",
					     0.0, 1.0, 0.0, 0.0,
					     0.0, 1.0, 0.0, 0.0,
					     0.0, 1.0);
      vacSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					   diffCoeffs, specCoeffs, GHz, GHz, GHz);

      //For the the interface of the Al and Cu    
      copTopSurfProp = new G4CMPSurfaceProperty("CopAlSurf",
						0.0, 1.0, 0.0, 0.0,
						0.0, 0.0, 0.0, 0.0,
						0.0, 1.0);
      copTopSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					      diffCoeffs, specCoeffs, GHz, GHz, GHz);

      //For the the interface of the Nb and Cu    
      copTopSurfProp2 = new G4CMPSurfaceProperty("CopNbSurf",
						 0.0, 1.0, 0.0, 0.0,
						 0.0, 0.0, 0.0, 0.0,
						 0.0, 1.0);
      copTopSurfProp2->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					       diffCoeffs, specCoeffs, GHz, GHz, GHz);

      
      //For the the interface of the Al and Si    
      topSurfProp = new G4CMPSurfaceProperty("AlSiSurf",
					     0.0, 1.0, 0.0, 0.0,
					     1.0, 0.0, 0.0, 0.0,
					     0.0, 1.0);
      topSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					   diffCoeffs, specCoeffs, GHz, GHz, GHz);

      //For the the interface of the Nb and Si
      topSurfProp2 = new G4CMPSurfaceProperty("NbSiSurf",
					      0.0, 1.0, 0.0, 0.0,
					      0.0, 0.0, 0.0, 0.0,
					      0.0, 1.0);
      topSurfProp2->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					    diffCoeffs, specCoeffs, GHz, GHz, GHz);
      
      //For the aluminum/niobium interface
      alNbSurfProp = new G4CMPSurfaceProperty("AlNbSurf",
					      0.0,1.0,0.0,0.0,
					      0.0,0.0,0.0,0.0,
					      0.0,0.0);
      alNbSurfProp->AddScatteringProperties(anhCutoff,reflCutoff,anhCoeffs,
					    diffCoeffs,specCoeffs,GHz,GHz,GHz);
      
      //For the interface of anything with the world
      wallSurfProp = new G4CMPSurfaceProperty("WallSurf",
					      0.0, 1.0, 0.0, 0.0,
					      0.0, 1.0, 0.0, 0.0,
					      0.0, 1.0);
      wallSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					    diffCoeffs, specCoeffs, GHz, GHz,GHz);

      
      
    }


    //
    // Separate surfaces for sensors vs. bare sidewall
    //
    new G4CMPLogicalBorderSurface("AlSi", SiPhys, aluminumTopPhysical,topSurfProp);
    new G4CMPLogicalBorderSurface("SiAl", aluminumTopPhysical, SiPhys, topSurfProp);
    new G4CMPLogicalBorderSurface("NbSi", SiPhys, niobiumTopPhysical,topSurfProp2);
    new G4CMPLogicalBorderSurface("SiNb", niobiumTopPhysical, SiPhys, topSurfProp2);

    new G4CMPLogicalBorderSurface("NbAl", niobiumTopPhysical,aluminumTopPhysical,alNbSurfProp);
    new G4CMPLogicalBorderSurface("AlNb", aluminumTopPhysical,niobiumTopPhysical,alNbSurfProp);


    
    new G4CMPLogicalBorderSurface("CuAl", copperOverTopPhysical, aluminumTopPhysical,copTopSurfProp);
    new G4CMPLogicalBorderSurface("AlCu", aluminumTopPhysical, copperOverTopPhysical,copTopSurfProp);
    new G4CMPLogicalBorderSurface("CuNb", copperOverTopPhysical, niobiumTopPhysical,copTopSurfProp2);
    new G4CMPLogicalBorderSurface("NbCu", niobiumTopPhysical,copperOverTopPhysical,copTopSurfProp2);
    
    new G4CMPLogicalBorderSurface("AlVac",aluminumTopPhysical,fWorldPhys,vacSurfProp);
    new G4CMPLogicalBorderSurface("VacAl",fWorldPhys,aluminumTopPhysical,vacSurfProp);
    new G4CMPLogicalBorderSurface("NbVac",niobiumTopPhysical,fWorldPhys,vacSurfProp);
    new G4CMPLogicalBorderSurface("VacNb",fWorldPhys,niobiumTopPhysical,vacSurfProp);

    
    new G4CMPLogicalBorderSurface("CuVac",copperOverTopPhysical,fWorldPhys,vacSurfProp);
    new G4CMPLogicalBorderSurface("VacCu",fWorldPhys,copperOverTopPhysical,vacSurfProp);
  
    //new G4CMPLogicalBorderSurface("VacAl",fWorldPhys,aluminumTopPhysical,vacSurfProp); //With this volume we try to reflect twice at each interface. BAD?
  
    new G4CMPLogicalBorderSurface("VacSi",SiPhys,fWorldPhys,wallSurfProp);
    new G4CMPLogicalBorderSurface("SiVac",fWorldPhys,SiPhys,wallSurfProp);

  
    //                                        
    // Visualization attributes
    //
    worldLogical->SetVisAttributes(G4VisAttributes::Invisible);
    G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    simpleBoxVisAtt->SetVisibility(true);
    siliconLogical->SetVisAttributes(simpleBoxVisAtt);
    aluminumLogical->SetVisAttributes(simpleBoxVisAtt);
    niobiumLogical->SetVisAttributes(simpleBoxVisAtt);
  }

  //-----------------------------------------------------------------------------------------------------------------------
  //Geometry 6:
  if( geometryID == 6 ){
    
    //                               
    // Silicon cylinder - this is the volume in which we will propagate phonons
    //  
    G4VSolid* siliconSolid = new G4Tubs("siliconSolid",0.*cm,3.81*cm,1.27*cm, 0.*deg, 360.*deg);
    G4LogicalVolume* siliconLogical = new G4LogicalVolume(siliconSolid,fSilicon,"siliconLogical");
    G4VPhysicalVolume* SiPhys = new G4PVPlacement(0,G4ThreeVector(),siliconLogical,"siliconPhysical", worldLogical,false,0);

  
    //
    // Niobium - large circle
    //
    G4VSolid* niobiumSolid = new G4Tubs("niobiumSolid",0.*cm,3.81*cm,4000000*nm, 0.*deg, 360.*deg);
    G4LogicalVolume* niobiumLogical = new G4LogicalVolume(niobiumSolid,fNiobium,"niobiumLogical");
    G4VPhysicalVolume* niobiumTopPhysical = new G4PVPlacement(0,G4ThreeVector(0.,0.,(1.27*cm + 4000000 * nm)), niobiumLogical,
							       "niobiumPhysical",worldLogical,false,0);

    
    //
    // Aluminum - smaller off-center circle
    //
    G4VSolid* aluminumSolid = new G4Tubs("aluminumSolid",0.*cm,1.0*mm,4000000*nm,0.*deg,360.*deg);
    G4LogicalVolume* aluminumLogical = new G4LogicalVolume(aluminumSolid,fAluminum,"aluminumLogical");
    G4VPhysicalVolume* aluminumTopPhysical = new G4PVPlacement(0,G4ThreeVector(1.0*cm,1.0*cm,0.0), aluminumLogical,
							       "aluminumPhysical",niobiumLogical,false,0);

    //
    // Niobium 2 - internal quarter-tubs with concavity
    //
    G4VSolid* niobiumSolid2 = new G4Tubs("niobiumSolid2",0.04*cm,0.05*cm,4000000*nm, 0.*deg, 90.*deg);
    G4LogicalVolume* niobiumLogical2 = new G4LogicalVolume(niobiumSolid2,fNiobium,"niobiumLogical2");
    G4VPhysicalVolume* niobiumTopPhysical2 = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.), niobiumLogical2,
							       "niobiumPhysical2",aluminumLogical,false,0);


    //
    // Niobium 3 - internal quarter-tubs with concavity
    //
    G4VSolid* niobiumSolid3 = new G4Tubs("niobiumSolid3",0.04*cm,0.05*cm,4000000*nm, 0.0*deg, 90.*deg);
    G4LogicalVolume* niobiumLogical3 = new G4LogicalVolume(niobiumSolid3,fNiobium,"niobiumLogical3");
    G4VPhysicalVolume* niobiumTopPhysical3 = new G4PVPlacement(0,G4ThreeVector(-0.3*mm,-0.3*mm,0.), niobiumLogical3,
							       "niobiumPhysical3",aluminumLogical,false,0);

  
    //
    // Copper, set to the full circle and above the Nb/Al layers
    //
    G4VSolid* copperSolid = new G4Tubs("aluminiumSolid",0.*cm,3.81*cm,0.01*cm,0.*deg, 360.*deg);
    G4LogicalVolume* copperLogical = new G4LogicalVolume(copperSolid,fCopper,"copperLogical");
    G4VPhysicalVolume* copperOverTopPhysical = new G4PVPlacement(0,G4ThreeVector(0.,0.,(1.27*cm + 8000000 * nm + 0.01 * cm)), copperLogical,
								 "copperLogical",worldLogical,false,0);


  
    //
    //Silicon lattice information
    //

    // G4LatticeManager gives physics processes access to lattices by volume
    G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();
    G4LatticeLogical* SiLogical = LM->LoadLattice(fSilicon, "Si");
    G4LatticeLogical* AlLogical = LM->LoadLattice(fAluminum, "Al");
    G4LatticeLogical* CuLogical = LM->LoadLattice(fCopper, "Cu");
    G4LatticeLogical* NbLogical = LM->LoadLattice(fNiobium, "Nb");
    G4LatticeLogical* NbLogical2 = LM->LoadLattice(fNiobium, "Nb");
    G4LatticeLogical* NbLogical3 = LM->LoadLattice(fNiobium, "Nb");
    
    // G4LatticePhysical assigns G4LatticeLogical a physical orientation
    G4LatticePhysical* SiPhysical = new G4LatticePhysical(SiLogical);
    G4LatticePhysical* AlPhysical = new G4LatticePhysical(AlLogical);
    G4LatticePhysical* CuPhysical = new G4LatticePhysical(CuLogical);
    G4LatticePhysical* NbPhysical = new G4LatticePhysical(NbLogical);
    G4LatticePhysical* NbPhysical2 = new G4LatticePhysical(NbLogical2);
    G4LatticePhysical* NbPhysical3 = new G4LatticePhysical(NbLogical3);
    SiPhysical->SetMillerOrientation(1,0,0);
    AlPhysical->SetMillerOrientation(1,0,0);
    CuPhysical->SetMillerOrientation(1,0,0);
    NbPhysical->SetMillerOrientation(1,0,0);
    NbPhysical2->SetMillerOrientation(1,0,0);
    NbPhysical3->SetMillerOrientation(1,0,0);
    LM->RegisterLattice(SiPhys, SiPhysical);
    LM->RegisterLattice(aluminumTopPhysical, AlPhysical);
    LM->RegisterLattice(copperOverTopPhysical, CuPhysical);
    LM->RegisterLattice(niobiumTopPhysical,NbPhysical);
    LM->RegisterLattice(niobiumTopPhysical2,NbPhysical2);
    LM->RegisterLattice(niobiumTopPhysical3,NbPhysical3);

    G4cout << "Al lattice LDOS: " << AlPhysical->GetLDOS() << ", scat rate: " << AlPhysical->GetPolycrystalElasticScatterMFP() << G4endl;

    //
    // surface between Al and Si determines phonon reflection/absorption
    //
    if (!fConstructed) {
      const G4double GHz = 1e9 * hertz; 

      //the following coefficients and cutoff values are not well-motivated
      //the code below is used only to demonstrate how to set these values.
      const std::vector<G4double> anhCoeffs = {0, 0, 0, 0, 0, 1.51e-14};
      const std::vector<G4double> diffCoeffs = {};
      const std::vector<G4double> specCoeffs = {};
      
    
      const G4double anhCutoff = 520., reflCutoff = 350.;   // Units external
      
      //For the the interface of the Cu and the world
      vacSurfProp = new G4CMPSurfaceProperty("VacCopSurf",
					     0.0, 1.0, 0.0, 0.0,
					     0.0, 1.0, 0.0, 0.0,
					     0.0, 1.0);
      vacSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					   diffCoeffs, specCoeffs, GHz, GHz, GHz);

      //For the the interface of the Al and Cu    
      copTopSurfProp = new G4CMPSurfaceProperty("CopAlSurf",
						0.0, 1.0, 0.0, 0.0,
						0.0, 0.0, 0.0, 0.0,
						0.0, 1.0);
      copTopSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					      diffCoeffs, specCoeffs, GHz, GHz, GHz);

      //For the the interface of the Nb and Cu    
      copTopSurfProp2 = new G4CMPSurfaceProperty("CopNbSurf",
						 0.0, 1.0, 0.0, 0.0,
						 0.0, 0.0, 0.0, 0.0,
						 0.0, 1.0);
      copTopSurfProp2->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					       diffCoeffs, specCoeffs, GHz, GHz, GHz);

      
      //For the the interface of the Al and Si    
      topSurfProp = new G4CMPSurfaceProperty("AlSiSurf",
					     0.0, 1.0, 0.0, 0.0,
					     1.0, 0.0, 0.0, 0.0,
					     0.0, 1.0);
      topSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					   diffCoeffs, specCoeffs, GHz, GHz, GHz);

      //For the the interface of the Nb and Si
      topSurfProp2 = new G4CMPSurfaceProperty("NbSiSurf",
					      0.0, 1.0, 0.0, 0.0,
					      1.0, 0.0, 0.0, 0.0,
					      0.0, 1.0);
      topSurfProp2->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					    diffCoeffs, specCoeffs, GHz, GHz, GHz);
      
      //For the aluminum/niobium interface
      alNbSurfProp = new G4CMPSurfaceProperty("AlNbSurf",
					      0.0,1.0,0.0,0.0,
					      0.0,0.0,0.0,0.0,
					      0.0,0.0);
      alNbSurfProp->AddScatteringProperties(anhCutoff,reflCutoff,anhCoeffs,
					    diffCoeffs,specCoeffs,GHz,GHz,GHz);
      
      //For the interface of anything with the world
      wallSurfProp = new G4CMPSurfaceProperty("WallSurf",
					      0.0, 1.0, 0.0, 0.0,
					      0.0, 1.0, 0.0, 0.0,
					      0.0, 1.0);
      wallSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					    diffCoeffs, specCoeffs, GHz, GHz,GHz);

      
      
    }


    //
    // Separate surfaces for sensors vs. bare sidewall
    //
    new G4CMPLogicalBorderSurface("AlSi", SiPhys, aluminumTopPhysical,topSurfProp);
    new G4CMPLogicalBorderSurface("SiAl", aluminumTopPhysical, SiPhys, topSurfProp);
    new G4CMPLogicalBorderSurface("NbSi", SiPhys, niobiumTopPhysical,topSurfProp2);
    new G4CMPLogicalBorderSurface("SiNb", niobiumTopPhysical, SiPhys, topSurfProp2);

    new G4CMPLogicalBorderSurface("NbAl", niobiumTopPhysical,aluminumTopPhysical,alNbSurfProp);
    new G4CMPLogicalBorderSurface("AlNb", aluminumTopPhysical,niobiumTopPhysical,alNbSurfProp);
    new G4CMPLogicalBorderSurface("Nb2Al", niobiumTopPhysical2,aluminumTopPhysical,alNbSurfProp);
    new G4CMPLogicalBorderSurface("AlNb2", aluminumTopPhysical,niobiumTopPhysical2,alNbSurfProp);
    new G4CMPLogicalBorderSurface("Nb3Al", niobiumTopPhysical3,aluminumTopPhysical,alNbSurfProp);
    new G4CMPLogicalBorderSurface("AlNb3", aluminumTopPhysical,niobiumTopPhysical3,alNbSurfProp);

    

    
    new G4CMPLogicalBorderSurface("CuAl", copperOverTopPhysical, aluminumTopPhysical,copTopSurfProp);
    new G4CMPLogicalBorderSurface("AlCu", aluminumTopPhysical, copperOverTopPhysical,copTopSurfProp);
    new G4CMPLogicalBorderSurface("CuNb", copperOverTopPhysical, niobiumTopPhysical,copTopSurfProp2);
    new G4CMPLogicalBorderSurface("NbCu", niobiumTopPhysical,copperOverTopPhysical,copTopSurfProp2);
    new G4CMPLogicalBorderSurface("CuNb2", copperOverTopPhysical, niobiumTopPhysical2,copTopSurfProp2);
    new G4CMPLogicalBorderSurface("Nb2Cu", niobiumTopPhysical2,copperOverTopPhysical,copTopSurfProp2);
    new G4CMPLogicalBorderSurface("CuNb3", copperOverTopPhysical, niobiumTopPhysical3,copTopSurfProp2);
    new G4CMPLogicalBorderSurface("Nb3Cu", niobiumTopPhysical3,copperOverTopPhysical,copTopSurfProp2);


    
    new G4CMPLogicalBorderSurface("AlVac",aluminumTopPhysical,fWorldPhys,vacSurfProp);
    new G4CMPLogicalBorderSurface("VacAl",fWorldPhys,aluminumTopPhysical,vacSurfProp);
    new G4CMPLogicalBorderSurface("NbVac",niobiumTopPhysical,fWorldPhys,vacSurfProp);
    new G4CMPLogicalBorderSurface("VacNb",fWorldPhys,niobiumTopPhysical,vacSurfProp);

    
    new G4CMPLogicalBorderSurface("CuVac",copperOverTopPhysical,fWorldPhys,vacSurfProp);
    new G4CMPLogicalBorderSurface("VacCu",fWorldPhys,copperOverTopPhysical,vacSurfProp);
  
    //new G4CMPLogicalBorderSurface("VacAl",fWorldPhys,aluminumTopPhysical,vacSurfProp); //With this volume we try to reflect twice at each interface. BAD?
  
    new G4CMPLogicalBorderSurface("VacSi",SiPhys,fWorldPhys,wallSurfProp);
    new G4CMPLogicalBorderSurface("SiVac",fWorldPhys,SiPhys,wallSurfProp);

  
    //                                        
    // Visualization attributes
    //
    worldLogical->SetVisAttributes(G4VisAttributes::Invisible);
    G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    simpleBoxVisAtt->SetVisibility(true);
    siliconLogical->SetVisAttributes(simpleBoxVisAtt);
    aluminumLogical->SetVisAttributes(simpleBoxVisAtt);
    niobiumLogical->SetVisAttributes(simpleBoxVisAtt);
    niobiumLogical2->SetVisAttributes(simpleBoxVisAtt);
    niobiumLogical3->SetVisAttributes(simpleBoxVisAtt);
  }
  
  else if( geometryID == 7 ){


    //Start by defining interface properties (since these are needed by classes we instantiate objects of here)
    if (!fConstructed) {
      const G4double GHz = 1e9 * hertz; 
      
      //the following coefficients and cutoff values are not well-motivated
      //the code below is used only to demonstrate how to set these values.
      const std::vector<G4double> anhCoeffs = {0, 0, 0, 0, 0, 1.51e-14};
      const std::vector<G4double> diffCoeffs = {};
      const std::vector<G4double> specCoeffs = {};
      
      
      const G4double anhCutoff = 520., reflCutoff = 350.;   // Units external
      
      //For the the interface of the Si and Aluminum
      fSiAlInterface = new G4CMPSurfaceProperty("SiAlSurf",
						0.0, 1.0, 0.0, 0.0,
						0.0, 0.0, 0.0, 0.0,
						0.0, 1.0);
      fSiAlInterface->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					      diffCoeffs, specCoeffs, GHz, GHz, GHz);
      fBorderContainer.emplace("SiAl",fSiAlInterface);

      //For the the interface of the Si and the world
      fSiVacInterface = new G4CMPSurfaceProperty("SiVacSurf",
						 0.0, 1.0, 0.0, 0.0,
						 0.0, 1.0, 0.0, 0.0,
						 0.0, 1.0);
      fSiVacInterface->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					       diffCoeffs, specCoeffs, GHz, GHz, GHz);
      fBorderContainer.emplace("SiVac",fSiVacInterface);

      //For the the interface of the Si and the Cu
      fSiCuInterface = new G4CMPSurfaceProperty("SiCuSurf",
						0.0, 1.0, 0.0, 0.0,
						0.0, 0.0, 0.0, 0.0,
						0.0, 1.0);
      fSiCuInterface->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					      diffCoeffs, specCoeffs, GHz, GHz, GHz);
      fBorderContainer.emplace("SiCu",fSiCuInterface);      

      //For the the interface of the Cu and the Vac
      fCuVacInterface = new G4CMPSurfaceProperty("CuVacSurf",
						 0.0, 1.0, 0.0, 0.0,
						 1.0, 1.0, 0.0, 0.0,
						 0.0, 1.0);
      fCuVacInterface->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					      diffCoeffs, specCoeffs, GHz, GHz, GHz);
      fBorderContainer.emplace("CuVac",fCuVacInterface);      

      
      //For the the interface of the Al and world
      fAlVacInterface = new G4CMPSurfaceProperty("AlVacSurf",
						 0.0, 1.0, 0.0, 0.0,
						 0.0, 1.0, 0.0, 0.0,
						 0.0, 1.0);
      fAlVacInterface->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					       diffCoeffs, specCoeffs, GHz, GHz, GHz);
      fBorderContainer.emplace("AlVac",fAlVacInterface);      

      //For the the interface of the Al and Al
      fAlAlInterface = new G4CMPSurfaceProperty("AlAlSurf",
						0.0, 1.0, 0.0, 0.0,
						0.0, 1.0, 0.0, 0.0,
						0.0, 0.0);
      fAlAlInterface->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					      diffCoeffs, specCoeffs, GHz, GHz, GHz);
      fBorderContainer.emplace("AlAl",fAlAlInterface);

      //For the the interface of the Vac and Vac
      fVacVacInterface = new G4CMPSurfaceProperty("VacVacSurf",
						  0.0, 1.0, 0.0, 0.0,
						  0.0, 1.0, 0.0, 0.0,
						  0.0, 1.0);
      fVacVacInterface->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
						diffCoeffs, specCoeffs, GHz, GHz, GHz);
      fBorderContainer.emplace("VacVac",fVacVacInterface);
      

    }

    //Also need to define logical lattices *here* now, since the logical lattice container needs to be passed into some classes
    // G4LatticeManager gives physics processes access to lattices by volume
    G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();
    G4LatticeLogical* log_siliconLattice = LM->LoadLattice(fSilicon, "Si");
    G4LatticeLogical* log_aluminumLattice = LM->LoadLattice(fAluminum, "Al");
    G4LatticeLogical* log_copperLattice = LM->LoadLattice(fCopper, "Cu");
    fLogicalLatticeContainer.emplace("Silicon",log_siliconLattice);
    fLogicalLatticeContainer.emplace("Aluminum",log_aluminumLattice);
    fLogicalLatticeContainer.emplace("Copper",log_copperLattice);

    

    
    //---------------------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------------
    // Now we start constructing the various components and their interfaces  
    bool checkOverlaps = true;


    //-------------------------------------------------------------------------------------------------------------------
    //First, set up the qubit chip substrate. 
    G4Box * solid_siliconChip = new G4Box("QubitChip_solid",
					  0.5*dp_siliconChipDimX,
					  0.5*dp_siliconChipDimY,
					  0.5*dp_siliconChipDimZ);
  
    //Now attribute a physical material to the chip
    G4LogicalVolume * log_siliconChip = new G4LogicalVolume(solid_siliconChip,
							    fSilicon,
							    "SiliconChip_log");
    
    //Now, create a physical volume and G4PVPlacement for storing as the final output
    G4ThreeVector siliconChipTranslate(0,0,0.5*(dp_housingDimZ - dp_siliconChipDimZ) + dp_eps); 
    G4VPhysicalVolume * phys_siliconChip = new G4PVPlacement(0,
							     siliconChipTranslate,
							     log_siliconChip,
							     "SiliconChip", 
							     worldLogical,
							     false,
							     0,
							     checkOverlaps);
    
    G4VisAttributes* siliconChipVisAtt= new G4VisAttributes(G4Colour(0.5,0.5,0.5));
    siliconChipVisAtt->SetVisibility(true);
    log_siliconChip->SetVisAttributes(siliconChipVisAtt);



    // G4LatticePhysical assigns G4LatticeLogical a physical orientation
    G4LatticePhysical* phys_siliconLattice = new G4LatticePhysical(log_siliconLattice);
    phys_siliconLattice->SetMillerOrientation(1,0,0); 
    LM->RegisterLattice(phys_siliconChip,phys_siliconLattice);

    //Set up border surfaces
    G4CMPLogicalBorderSurface * border_siliconChip_world = new G4CMPLogicalBorderSurface("border_siliconChip_world", phys_siliconChip, fWorldPhys, fSiVacInterface);

    



    //-------------------------------------------------------------------------------------------------------------------
    //If desired, set up the copper qubit housing
    if( dp_useQubitHousing ){
      
      
      QuasiparticleQubitHousing * qubitHousing = new QuasiparticleQubitHousing(0,
									       G4ThreeVector(0,0,0),
									       "QubitHousing",
									       worldLogical,
									       false,
									       0,
									       checkOverlaps);
      G4LogicalVolume * log_qubitHousing = qubitHousing->GetLogicalVolume();
      G4VPhysicalVolume * phys_qubitHousing = qubitHousing->GetPhysicalVolume();

      //Set up lattice properties for copper housing
      G4LatticePhysical* phys_copperLattice = new G4LatticePhysical(log_copperLattice);
      phys_copperLattice->SetMillerOrientation(1,0,0); //Silly, but let's do it for kicks
      LM->RegisterLattice(phys_qubitHousing,phys_copperLattice);


      
      //Set up the logical border surface
      G4CMPLogicalBorderSurface * border_siliconChip_qubitHousing = new G4CMPLogicalBorderSurface("border_siliconChip_qubitHousing", phys_siliconChip, phys_qubitHousing, fSiCuInterface);
      G4CMPLogicalBorderSurface * border_qubitHousing_siliconChip = new G4CMPLogicalBorderSurface("border_qubitHousing_siliconChip", phys_qubitHousing, phys_siliconChip, fSiCuInterface);
      G4CMPLogicalBorderSurface * border_qubitHousing_world = new G4CMPLogicalBorderSurface("border_qubitHousing_world", phys_qubitHousing, fWorldPhys,fCuVacInterface);

    }
    


    
   
    //-------------------------------------------------------------------------------------------------------------------
    //Now set up the ground plane, in which the transmission line, resonators, and qubits will be located.
    if( dp_useGroundPlane ){
    
    
      G4Box * solid_groundPlane = new G4Box("GroundPlane_solid",
					    0.5*dp_groundPlaneDimX,
					    0.5*dp_groundPlaneDimY,
					    0.5*dp_groundPlaneDimZ);
    
    
      //Now attribute a physical material to the chip
      G4LogicalVolume * log_groundPlane = new G4LogicalVolume(solid_groundPlane,
							      fNiobium,
							      "GroundPlane_log");
    
    
      //Now, create a physical volume and G4PVPlacement for storing as the final output
      G4ThreeVector groundPlaneTranslate(0,0,0.5*(dp_housingDimZ) + dp_eps + dp_groundPlaneDimZ*0.5);
      G4VPhysicalVolume * phys_groundPlane = new G4PVPlacement(0,
							       groundPlaneTranslate,
							       log_groundPlane,
							       "GroundPlane", 
							       worldLogical,
							       false,
							       0,
							       checkOverlaps);
      
      G4VisAttributes* groundPlaneVisAtt= new G4VisAttributes(G4Colour(0.0,1.0,1.0,0.5));
      groundPlaneVisAtt->SetVisibility(true);
      log_groundPlane->SetVisAttributes(groundPlaneVisAtt);

      G4LatticePhysical* phys_groundPlaneLattice = new G4LatticePhysical(log_aluminumLattice);
      phys_groundPlaneLattice->SetMillerOrientation(1,0,0); //Silly for a polycrystalline thing, but okay
      LM->RegisterLattice(phys_groundPlane,phys_groundPlaneLattice);

      
      
      //Set up the logical border surface
      G4CMPLogicalBorderSurface * border_siliconChip_groundPlane = new G4CMPLogicalBorderSurface("border_siliconChip_groundPlane", phys_siliconChip, phys_groundPlane, fSiAlInterface);
      G4CMPLogicalBorderSurface * border_groundPlane_siliconChip = new G4CMPLogicalBorderSurface("border_siliconChip_groundPlane", phys_groundPlane, phys_siliconChip, fSiAlInterface);
      G4CMPLogicalBorderSurface * border_world_groundPlane = new G4CMPLogicalBorderSurface("border_world_groundPlane", fWorldPhys, phys_groundPlane, fAlVacInterface);
      G4CMPLogicalBorderSurface * border_groundPlane_world = new G4CMPLogicalBorderSurface("border_groundPlane_world", phys_groundPlane, fWorldPhys, fAlVacInterface);


    
    
      //-------------------------------------------------------------------------------------------------------------------
      //Now set up the transmission line
      if( dp_useTransmissionLine ){
	
	G4ThreeVector transmissionLineTranslate(0,0,0.0);//Since it's within the ground plane exactly; 0.5*(dp_housingDimZ) + dp_eps + dp_groundPlaneDimZ*0.5 ); 
	QuasiparticleTransmissionLine * tLine = new QuasiparticleTransmissionLine(0,
										  transmissionLineTranslate,
										  "TransmissionLine",
										  log_groundPlane,
										  false,
										  0,
										  LM,
										  fLogicalLatticeContainer,
										  fBorderContainer,
										  checkOverlaps);
	G4LogicalVolume * log_tLine = tLine->GetLogicalVolume();
	G4VPhysicalVolume * phys_tLine = tLine->GetPhysicalVolume();



	//Now, if we're using the chip and ground plane AND the transmission line
	//This gets a bit hairy, since the transmission line is composite of both Nb and vacuum.
	//So we'll access the list of physical objects present in it and link those one-by-one to the
	//silicon chip and to the world
	for( int iSubVol = 0; iSubVol < tLine->GetListOfAllFundamentalSubVolumes().size(); ++iSubVol){
	  std::cout << "TLine sub volume names (to be used for boundaries): " << std::get<1>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]) << " with material " << std::get<0>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]) << std::endl;

	  //Set the si/vac and si/Al interfaces to be symmetric in both dimensions
	  std::string tempName1a = "border_siliconChip_" + std::get<1>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]);
	  std::string tempName1b = "border_" + std::get<1>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]) + "_siliconChip";	  
	  if( std::get<0>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Vacuum") != std::string::npos ){
	    G4CMPLogicalBorderSurface * border1_siliconChip_transmissionLineEmpty = new G4CMPLogicalBorderSurface(tempName1a, phys_siliconChip, std::get<2>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fSiVacInterface);
	    G4CMPLogicalBorderSurface * border2_siliconChip_transmissionLineEmpty = new G4CMPLogicalBorderSurface(tempName1b, std::get<2>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), phys_siliconChip, fSiVacInterface);
	  }
	  if( std::get<0>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Aluminum") != std::string::npos ){
	    G4CMPLogicalBorderSurface * border1_siliconChip_transmissionLineConductor = new G4CMPLogicalBorderSurface(tempName1a, phys_siliconChip, std::get<2>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fSiAlInterface);
	    G4CMPLogicalBorderSurface * border2_siliconChip_transmissionLineConductor = new G4CMPLogicalBorderSurface(tempName1b, std::get<2>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), phys_siliconChip, fSiAlInterface);
	  }


	  //Set the world/vac and world/Al interfaces to be symmetric in both dimensions
	  std::string tempName2a = "border_world_" + std::get<1>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]);
	  std::string tempName2b = "border_" + std::get<1>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]) + "_world";	  
	  if( std::get<0>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Vacuum") != std::string::npos ){
	    G4CMPLogicalBorderSurface * border1_world_transmissionLineEmpty = new G4CMPLogicalBorderSurface(tempName2a, fWorldPhys, std::get<2>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fVacVacInterface);
	    G4CMPLogicalBorderSurface * border2_world_transmissionLineEmpty = new G4CMPLogicalBorderSurface(tempName2b, std::get<2>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fWorldPhys, fVacVacInterface);
	  }
	  if( std::get<0>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Aluminum") != std::string::npos ){
	    G4CMPLogicalBorderSurface * border1_world_transmissionLineConductor = new G4CMPLogicalBorderSurface(tempName2a, fWorldPhys, std::get<2>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fAlVacInterface);
	    G4CMPLogicalBorderSurface * border2_world_transmissionLineConductor = new G4CMPLogicalBorderSurface(tempName2b, std::get<2>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fWorldPhys, fAlVacInterface);
	  }


	  //Space here for linking to GROUND PLANE. Only need to link vacuum ones since they're fully enclosing TL
	  std::string tempName3a = "border_groundPlane_" + std::get<1>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]);
	  std::string tempName3b = "border_" + std::get<1>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]) + "_groundPlane";	  
	  if( std::get<0>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Vacuum") != std::string::npos ){
	    G4CMPLogicalBorderSurface * border1_groundPlane_transmissionLineEmpty = new G4CMPLogicalBorderSurface(tempName3a, phys_groundPlane, std::get<2>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fAlVacInterface);
	    G4CMPLogicalBorderSurface * border2_groundPlane_transmissionLineEmpty = new G4CMPLogicalBorderSurface(tempName3b, std::get<2>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), phys_groundPlane, fAlVacInterface);
	  }
	  
	}
      }
    }//REMOVE
    /*

      //-------------------------------------------------------------------------------------------------------------------
      //Now set up a set of 6 resonator assemblies
      if( dp_useResonatorAssembly ){
	int nR = 6;
	for( int iR = 0; iR < nR; ++iR ){
      
	  //First, get the translation vector for the resonator assembly
	  //For the top three, don't do a rotation. For the bottom three, do
	  G4ThreeVector resonatorAssemblyTranslate(0,0,0);
	  G4RotationMatrix * rotAssembly = 0;
	  if( iR <= 2 ){
	    resonatorAssemblyTranslate = G4ThreeVector(dp_resonatorLateralSpacing*(iR-1)+dp_centralResonatorOffsetX,
						       0.5 * dp_resonatorAssemblyBaseNbDimY + 0.5 * dp_transmissionLineCavityFullWidth,
						       0.0);
	    rotAssembly = 0;
	  }
	  else{
	    resonatorAssemblyTranslate = G4ThreeVector(dp_resonatorLateralSpacing*(iR-4)-dp_centralResonatorOffsetX, //Negative offset because qubit is mirrored on underside
						       -1*(0.5 * dp_resonatorAssemblyBaseNbDimY + 0.5 * dp_transmissionLineCavityFullWidth),
						       0.0);
	    rotAssembly = new G4RotationMatrix();
	    rotAssembly->rotateZ(180*deg);
	  }
	
	  char name[400];
	  sprintf(name,"ResonatorAssembly_%d",iR);
	  G4String resonatorAssemblyName(name);
	  RISQTutorialResonatorAssembly * resonatorAssembly = new RISQTutorialResonatorAssembly(rotAssembly,
												resonatorAssemblyTranslate,
												resonatorAssemblyName,
												log_groundPlane,
												false,
												0,
												checkOverlaps);
	  G4LogicalVolume * log_resonatorAssembly = resonatorAssembly->GetLogicalVolume();
	  G4VPhysicalVolume * phys_resonatorAssembly = resonatorAssembly->GetPhysicalVolume();
	
	
	  //Do the logical border creation now
	  for( int iSubVol = 0; iSubVol < resonatorAssembly->GetListOfAllFundamentalSubVolumes().size(); ++iSubVol){
	    std::cout << "TLine sub volume names (to be used for boundaries): " << std::get<1>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]) << " with material " << std::get<0>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]) << std::endl;
	  
	    std::string tempName = "border_siliconChip_" + std::get<1>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]);	  
	    if( std::get<0>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Vacuum") != std::string::npos ){
	      G4CMPLogicalBorderSurface * border_siliconChip_resonatorAssemblyEmpty = new G4CMPLogicalBorderSurface(tempName, phys_siliconChip, std::get<2>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]), fSiVacuumInterface);
	    }
	    if( std::get<0>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Niobium") != std::string::npos ){
	      G4CMPLogicalBorderSurface * border_siliconChip_resonatorAssemblyConductor = new G4CMPLogicalBorderSurface(tempName, phys_siliconChip, std::get<2>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]), fSiNbInterface);
	    }
	  }
	}
      }
    
    
    
      //-------------------------------------------------------------------------------------------------------------------
      // Flux lines
      if( dp_useFluxLines ){
      
      
	//--------------------
	G4ThreeVector topStraightFluxLineTranslate(dp_topCenterFluxLineOffsetX,dp_topCenterFluxLineOffsetY,0);
	RISQTutorialStraightFluxLine * topStraightFLine = new RISQTutorialStraightFluxLine(0,
											   topStraightFluxLineTranslate,
											   "TopStraightFluxLine",
											   log_groundPlane,
											   false,
											   0,
											   checkOverlaps);
	G4LogicalVolume * log_topStraightFline = topStraightFLine->GetLogicalVolume();
	G4VPhysicalVolume * phys_topStraightFline = topStraightFLine->GetPhysicalVolume();
      
	//Do the logical border creation now
	for( int iSubVol = 0; iSubVol < topStraightFLine->GetListOfAllFundamentalSubVolumes().size(); ++iSubVol){
	  std::cout << "TLine sub volume names (to be used for boundaries): " << std::get<1>(topStraightFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]) << " with material " << std::get<0>(topStraightFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]) << std::endl;
	
	  std::string tempName = "border_siliconChip_" + std::get<1>(topStraightFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]);	  
	  if( std::get<0>(topStraightFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Vacuum") != std::string::npos ){
	    G4CMPLogicalBorderSurface * border_siliconChip_topStraightFluxLineEmpty = new G4CMPLogicalBorderSurface(tempName, phys_siliconChip, std::get<2>(topStraightFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fSiVacuumInterface);
	  }
	  if( std::get<0>(topStraightFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Niobium") != std::string::npos ){
	    G4CMPLogicalBorderSurface * border_siliconChip_topStraightFluxLineConductor = new G4CMPLogicalBorderSurface(tempName, phys_siliconChip, std::get<2>(topStraightFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fSiNbInterface);
	  }
	}




	


	//--------------------
	G4ThreeVector bottomStraightFluxLineTranslate(dp_topCenterFluxLineOffsetX,-1*dp_topCenterFluxLineOffsetY,0);
	G4RotationMatrix * rotBottomCenter = new G4RotationMatrix();
	rotBottomCenter->rotateZ(180.*deg);
	RISQTutorialStraightFluxLine * bottomStraightFLine = new RISQTutorialStraightFluxLine(rotBottomCenter,
											      bottomStraightFluxLineTranslate,
											      "BottomStraightFluxLine",
											      log_groundPlane,
											      false,
											      0,
											      checkOverlaps);
	G4LogicalVolume * log_bottomStraightFline = bottomStraightFLine->GetLogicalVolume();
	G4VPhysicalVolume * phys_bottomStraightFline = bottomStraightFLine->GetPhysicalVolume();

	//Do the logical border creation now
	for( int iSubVol = 0; iSubVol < bottomStraightFLine->GetListOfAllFundamentalSubVolumes().size(); ++iSubVol){
	  std::cout << "TLine sub volume names (to be used for boundaries): " << std::get<1>(bottomStraightFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]) << " with material " << std::get<0>(bottomStraightFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]) << std::endl;

	  std::string tempName = "border_siliconChip_" + std::get<1>(bottomStraightFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]);	  
	  if( std::get<0>(bottomStraightFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Vacuum") != std::string::npos ){
	    G4CMPLogicalBorderSurface * border_siliconChip_bottomStraightFluxLineEmpty = new G4CMPLogicalBorderSurface(tempName, phys_siliconChip, std::get<2>(bottomStraightFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fSiVacuumInterface);
	  }
	  if( std::get<0>(bottomStraightFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Niobium") != std::string::npos ){
	    G4CMPLogicalBorderSurface * border_siliconChip_bottomStraightFluxLineConductor = new G4CMPLogicalBorderSurface(tempName, phys_siliconChip, std::get<2>(bottomStraightFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fSiNbInterface);
	  }

	}

	

       

	//--------------------
	//Corner flux line 1
	G4ThreeVector topLeftCornerFluxLineTranslate(dp_topLeftFluxLineOffsetX,dp_topLeftFluxLineOffsetY,0);
	G4RotationMatrix * rotTopLeftCenter = new G4RotationMatrix();
	rotTopLeftCenter->rotateZ(0.*deg);
	RISQTutorialCornerFluxLine * topLeftCornerFLine = new RISQTutorialCornerFluxLine(rotTopLeftCenter,
											 topLeftCornerFluxLineTranslate,
											 "GroundPlane_TopLeftCornerFluxLine",
											 log_groundPlane,
											 false,
											 0,
											 checkOverlaps);


	
	G4LogicalVolume * log_topLeftCornerFline = topLeftCornerFLine->GetLogicalVolume();
	G4VPhysicalVolume * phys_topLeftCornerFline = topLeftCornerFLine->GetPhysicalVolume();

	
	//Do the logical border creation now
	for( int iSubVol = 0; iSubVol < topLeftCornerFLine->GetListOfAllFundamentalSubVolumes().size(); ++iSubVol){
	  std::cout << "TLine sub volume names (to be used for boundaries): " << std::get<1>(topLeftCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]) << " with material " << std::get<0>(topLeftCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]) << std::endl;

	  std::string tempName = "border_siliconChip_" + std::get<1>(topLeftCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]);	  
	  if( std::get<0>(topLeftCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Vacuum") != std::string::npos ){
	    G4CMPLogicalBorderSurface * border_siliconChip_topLeftCornerFluxLineEmpty = new G4CMPLogicalBorderSurface(tempName, phys_siliconChip, std::get<2>(topLeftCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fSiVacuumInterface);
	  }
	  if( std::get<0>(topLeftCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Niobium") != std::string::npos ){
	    G4CMPLogicalBorderSurface * border_siliconChip_topLeftCornerFluxLineConductor = new G4CMPLogicalBorderSurface(tempName, phys_siliconChip, std::get<2>(topLeftCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fSiNbInterface);
	  }	  
	}






	

	//--------------------
	//Corner flux line 2
	G4ThreeVector topRightCornerFluxLineTranslate(-1*dp_topLeftFluxLineOffsetX,dp_topLeftFluxLineOffsetY,0);
	G4RotationMatrix * rotTopRightCenter = new G4RotationMatrix();
	rotTopRightCenter->rotateY(180.*deg);
	RISQTutorialCornerFluxLine * topRightCornerFLine = new RISQTutorialCornerFluxLine(rotTopRightCenter,
											  topRightCornerFluxLineTranslate,
											  "GroundPlane_TopRightCornerFluxLine",
											  log_groundPlane,
											  false,
											  0,
											  checkOverlaps);
	G4LogicalVolume * log_topRightCornerFline = topRightCornerFLine->GetLogicalVolume();
	G4VPhysicalVolume * phys_topRightCornerFline = topRightCornerFLine->GetPhysicalVolume();

	//Do the logical border creation now
	for( int iSubVol = 0; iSubVol < topRightCornerFLine->GetListOfAllFundamentalSubVolumes().size(); ++iSubVol){
	  std::cout << "TLine sub volume names (to be used for boundaries): " << std::get<1>(topRightCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]) << " with material " << std::get<0>(topRightCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]) << std::endl;

	  std::string tempName = "border_siliconChip_" + std::get<1>(topRightCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]);	  
	  if( std::get<0>(topRightCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Vacuum") != std::string::npos ){
	    G4CMPLogicalBorderSurface * border_siliconChip_topRightCornerFluxLineEmpty = new G4CMPLogicalBorderSurface(tempName, phys_siliconChip, std::get<2>(topRightCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fSiVacuumInterface);
	  }
	  if( std::get<0>(topRightCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Niobium") != std::string::npos ){
	    G4CMPLogicalBorderSurface * border_siliconChip_topRightCornerFluxLineConductor = new G4CMPLogicalBorderSurface(tempName, phys_siliconChip, std::get<2>(topRightCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fSiNbInterface);
	  }

	}

    


	//--------------------
	//Corner flux line 3
	G4ThreeVector bottomLeftCornerFluxLineTranslate(dp_topLeftFluxLineOffsetX,-1*dp_topLeftFluxLineOffsetY,0);
	G4RotationMatrix * rotBottomLeftCenter = new G4RotationMatrix();
	rotBottomLeftCenter->rotateX(180.*deg);
	RISQTutorialCornerFluxLine * bottomLeftCornerFLine = new RISQTutorialCornerFluxLine(rotBottomLeftCenter,
											    bottomLeftCornerFluxLineTranslate,
											    "BottomLeftCornerFluxLine",
											    log_groundPlane,
											    false,
											    0,
											    checkOverlaps);
	G4LogicalVolume * log_bottomLeftCornerFline = bottomLeftCornerFLine->GetLogicalVolume();
	G4VPhysicalVolume * phys_bottomLeftCornerFline = bottomLeftCornerFLine->GetPhysicalVolume();

	//Do the logical border creation now
	for( int iSubVol = 0; iSubVol < bottomLeftCornerFLine->GetListOfAllFundamentalSubVolumes().size(); ++iSubVol){
	  std::cout << "TLine sub volume names (to be used for boundaries): " << std::get<1>(bottomLeftCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]) << " with material " << std::get<0>(bottomLeftCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]) << std::endl;

	  
	  std::string tempName = "border_siliconChip_" + std::get<1>(bottomLeftCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]);	  
	  if( std::get<0>(bottomLeftCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Vacuum") != std::string::npos ){
	    G4CMPLogicalBorderSurface * border_siliconChip_bottomLeftCornerFluxLineEmpty = new G4CMPLogicalBorderSurface(tempName, phys_siliconChip, std::get<2>(bottomLeftCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fSiVacuumInterface);
	  }
	  if( std::get<0>(bottomLeftCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Niobium") != std::string::npos ){
	    G4CMPLogicalBorderSurface * border_siliconChip_bottomLeftCornerFluxLineConductor = new G4CMPLogicalBorderSurface(tempName, phys_siliconChip, std::get<2>(bottomLeftCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fSiNbInterface);
	  }

	  
	}

    
    
	//--------------------
	//Corner flux line 4
	G4ThreeVector bottomRightCornerFluxLineTranslate(-1*dp_topLeftFluxLineOffsetX,-1*dp_topLeftFluxLineOffsetY,0);
	G4RotationMatrix * rotBottomRightCenter = new G4RotationMatrix();
	rotBottomRightCenter->rotateX(180.*deg);
	rotBottomRightCenter->rotateY(180.*deg);
	RISQTutorialCornerFluxLine * bottomRightCornerFLine = new RISQTutorialCornerFluxLine(rotBottomRightCenter,
											     bottomRightCornerFluxLineTranslate,
											     "BottomRightCornerFluxLine",
											     log_groundPlane,
											     false,
											     0,
											     checkOverlaps);
	G4LogicalVolume * log_bottomRightCornerFline = bottomRightCornerFLine->GetLogicalVolume();
	G4VPhysicalVolume * phys_bottomRightCornerFline = bottomRightCornerFLine->GetPhysicalVolume();

	//Do the logical border creation now
	for( int iSubVol = 0; iSubVol < bottomRightCornerFLine->GetListOfAllFundamentalSubVolumes().size(); ++iSubVol){
	  std::cout << "TLine sub volume names (to be used for boundaries): " << std::get<1>(bottomRightCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]) << " with material " << std::get<0>(bottomRightCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]) << std::endl;

	  std::string tempName = "border_siliconChip_" + std::get<1>(bottomRightCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]);	  
	  if( std::get<0>(bottomRightCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Vacuum") != std::string::npos ){
	    G4CMPLogicalBorderSurface * border_siliconChip_bottomRightCornerFluxLineEmpty = new G4CMPLogicalBorderSurface(tempName, phys_siliconChip, std::get<2>(bottomRightCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fSiVacuumInterface);
	  }
	  if( std::get<0>(bottomRightCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Niobium") != std::string::npos ){
	    G4CMPLogicalBorderSurface * border_siliconChip_bottomRightCornerFluxLineConductor = new G4CMPLogicalBorderSurface(tempName, phys_siliconChip, std::get<2>(bottomRightCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fSiNbInterface);
	  }
	}
      }
    }



    
    */
  }


  
}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Attach material properties and electrode/sensor handler to surface

void QuasiparticleDetectorConstruction::
AttachPhononSensor(G4CMPSurfaceProperty *surfProp) {
  if (!surfProp) return;		// No surface, nothing to do

  // Specify properties of aluminum sensor, same on both detector faces
  // See G4CMPPhononElectrode.hh or README.md for property keys

  /*
  // Properties must be added to existing surface-property table
  auto sensorProp = surfProp->GetPhononMaterialPropertiesTablePointer();
  sensorProp->AddConstProperty("filmAbsorption", 0.20);    // True sensor area
  sensorProp->AddConstProperty("filmThickness", 600.*nm);
  sensorProp->AddConstProperty("gapEnergy", 173.715e-6*eV);
  sensorProp->AddConstProperty("lowQPLimit", 3.);
  sensorProp->AddConstProperty("phononLifetime", 242.*ps);
  sensorProp->AddConstProperty("phononLifetimeSlope", 0.29);
  sensorProp->AddConstProperty("vSound", 3.26*km/s);
  sensorProp->AddConstProperty("subgapAbsorption", 0.1);

  // Attach electrode object to handle KaplanQP interface
  surfProp->SetPhononElectrode(new G4CMPPhononElectrode);
  */
}


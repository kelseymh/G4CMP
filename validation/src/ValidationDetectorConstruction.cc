/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file exoticphysics/phonon/src/ValidationDetectorConstruction.cc \brief
/// Implementation of the ValidationDetectorConstruction class
//
// $Id: a2016d29cc7d1e75482bfc623a533d20b60390da $
//
// 20140321  Drop passing placement transform to G4LatticePhysical
// 20211207  Replace G4Logical*Surface with G4CMP-specific versions.
// 20220809  [ For M. Hui ] -- Add frequency dependent surface properties.
// 20221006  Remove unused features; add phonon sensor pad with use of
//		G4CMPPhononElectrode to demonstrate KaplanQP.

#include "ValidationDetectorConstruction.hh"
#include "ValidationDetectorParameters.hh"
#include "ValidationSensitivity.hh"
#include "ValidationQubitHousing.hh"
#include "ValidationTransmissionLine.hh"
#include "ValidationResonatorAssembly.hh"
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
#include "ValidationConfigManager.hh"

using namespace ValidationDetectorParameters;


ValidationDetectorConstruction::ValidationDetectorConstruction()
  : fLiquidHelium(0), fSilicon(0), fAluminum(0), fTungsten(0), fCopper(0),
    fNiobium(0), fWorldPhys(0), fVacSurfProp(0), fSiGeSurfProp(0),
    fGeAl1SurfProp(0), fAl1Al2FromAboveSurfProp(0),fAl1Al2FromBelowSurfProp(0),
    fAl2NbSurfProp(0), fAl2Al3SurfProp(0), fAl3NbSurfProp(0),
    electrodeSensitivity(0), fConstructed(false) {;}


ValidationDetectorConstruction::~ValidationDetectorConstruction() {
  delete fVacSurfProp;
  delete fSiGeSurfProp;
  delete fGeAl1SurfProp;
  delete fAl1Al2FromAboveSurfProp;
  delete fAl1Al2FromBelowSurfProp;
  delete fAl2NbSurfProp;
  delete fAl2Al3SurfProp;
  delete fAl3NbSurfProp;
}

G4VPhysicalVolume* ValidationDetectorConstruction::Construct() {
  G4cout << "fConstructed first status: " << fConstructed << G4endl;
  
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

void ValidationDetectorConstruction::DefineMaterials() { 
  G4NistManager* nistManager = G4NistManager::Instance();

  fLiquidHelium = nistManager->FindOrBuildMaterial("G4_AIR"); // to be corrected
  fSilicon = nistManager->FindOrBuildMaterial("G4_Si");
  fGermanium = nistManager->FindOrBuildMaterial("G4_Ge");
  fAluminum = nistManager->FindOrBuildMaterial("G4_Al");
  fCopper = nistManager->FindOrBuildMaterial("G4_Cu");
  fTungsten = nistManager->FindOrBuildMaterial("G4_W");
  fNiobium = nistManager->FindOrBuildMaterial("G4_Nb");
}


void ValidationDetectorConstruction::SetupGeometry() {
  //Define the config manager
  int geometryID = ValidationConfigManager::Instance()->GetGeometryID();
  G4cout << "Starting with validation geometryID: " << geometryID << G4endl;
  G4cout << "Testing config manager call 2: "
         << ValidationConfigManager::GetGeometryID() << G4endl;
  
  //     
  // World
  //
  G4VSolid* worldSolid = new G4Box("World",16.*cm,16.*cm,16.*cm);
  G4LogicalVolume* worldLogical =
    new G4LogicalVolume(worldSolid,fLiquidHelium,"World");
  worldLogical->SetUserLimits(new G4UserLimits(100*mm, DBL_MAX, DBL_MAX, 0, 0));
  fWorldPhys = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"World",0,
                                 false,0);


  //For this validation script, we allow toggling between a few geometries.
  //1. Geometry 1: A multipurpose geometry for validating that SC/QP Tracking
  //   is working as intended
  //2. Geometry 2: A "full qubit chip" geometry for testing a comprehensive
  //   set of the physics processes together
  //3. Geometry X: You may add geometries to this if you want to add validation
  //   tests for those!
  
  //Geometry 1: Multipurpose phonon and QP tracking validation geometry
  if (geometryID == 1) {

    //Define volumes first
    
    //First, we define a non-superconductor cylinder: Si
    G4VSolid* siliconSolid = new G4Tubs("siliconSolid",0.*cm,1.0*cm,
                                        dp_siThickness/2.0, 0.*deg, 360.*deg);
    G4LogicalVolume* siliconLogical = new G4LogicalVolume(siliconSolid,fSilicon,
                                                          "siliconLogical");
    G4VPhysicalVolume* siliconPhysical
      = new G4PVPlacement(0,G4ThreeVector(),siliconLogical,"siliconPhysical",
                          worldLogical,false,0);
    
    //Next, we define another non-superconductor cylinder: Ge. Place it just
    //below the Si and in contact
    G4VSolid* germaniumSolid = new G4Tubs("germaniumSolid",0.*cm,1.0*cm,
                                          dp_geThickness/2.0, 0.*deg, 360.*deg);
    G4LogicalVolume* germaniumLogical
      = new G4LogicalVolume(germaniumSolid,fGermanium,"germaniumLogical");
    G4VPhysicalVolume* germaniumPhysical
      = new G4PVPlacement(0,
                          G4ThreeVector(0.0,0.0,(-dp_siThickness/2.0-dp_geThickness/2.0)),
                          germaniumLogical,"germaniumPhysical", worldLogical,
                          false,0);
    
    //Next, we define a superconductor cylinder: Al1. Place it just below the
    //Ge and in contact. Needs to be thin enough so that quasielastic
    //scattering tests don't take forever, with one of the top/bottom
    //boundaries being fully absorptive
    G4VSolid* aluminum1Solid
      = new G4Tubs("aluminum1Solid",0.*cm,1.0*cm,dp_aluminum1Thickness/2.0,
                   0.*deg,360.*deg);
    G4LogicalVolume* aluminum1Logical
      = new G4LogicalVolume(aluminum1Solid,fAluminum,"aluminum1Logical");
    G4VPhysicalVolume* aluminum1Physical
      = new G4PVPlacement(0,
                          G4ThreeVector(0.0,0.0,(-dp_siThickness/2.0 - dp_geThickness - dp_aluminum1Thickness/2.0)),
                          aluminum1Logical,"aluminum1Physical", worldLogical,
                          false,0);
    
    //Next, we define another superconductor, but make it a long G4Box: Al2.
    //This long strip tests diffusion physics. Place it just below Al1 and in
    //contact. Making it the same thickness as Al1 arbitrarily.
    G4VSolid* aluminum2Solid
      = new G4Box("aluminum2Solid",0.1*cm,1.0*cm,dp_aluminum2Thickness/2.0);
    G4LogicalVolume* aluminum2Logical
      = new G4LogicalVolume(aluminum2Solid,fAluminum,"aluminum2Logical");
    G4VPhysicalVolume* aluminum2Physical
      = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,(-dp_siThickness/2.0 - dp_geThickness - dp_aluminum1Thickness - dp_aluminum2Thickness/2.0)),
                          aluminum2Logical,"aluminum2Physical", worldLogical,
                          false,0);

    //We define two Nb blocks to flank aluminum2
    G4VSolid* niobiumAForAluminum2Solid
      = new G4Box("niobiumAForAluminum2Solid",0.1*cm,0.1*cm,
                  dp_aluminum2Thickness/2.0);
    G4LogicalVolume* niobiumAForAluminum2Logical
      = new G4LogicalVolume(niobiumAForAluminum2Solid,fNiobium,
                            "niobiumAForAluminum2Logical");
    G4VPhysicalVolume* niobiumAForAluminum2Physical
      = new G4PVPlacement(0,G4ThreeVector(0.0,-1.1*cm,(-dp_siThickness/2.0 - dp_geThickness - dp_aluminum1Thickness - dp_aluminum2Thickness/2.0)),
                          niobiumAForAluminum2Logical,
                          "niobiumAForAluminum2Physical", worldLogical,false,0);

    G4VSolid* niobiumBForAluminum2Solid
      = new G4Box("niobiumBForAluminum2Solid",0.1*cm,0.1*cm,
                  dp_aluminum2Thickness/2.0);
    G4LogicalVolume* niobiumBForAluminum2Logical
      = new G4LogicalVolume(niobiumBForAluminum2Solid,fNiobium,
                            "niobiumBForAluminum2Logical");
    G4VPhysicalVolume* niobiumBForAluminum2Physical
      = new G4PVPlacement(0,G4ThreeVector(0.0,1.1*cm,(-dp_siThickness/2.0 - dp_geThickness - dp_aluminum1Thickness - dp_aluminum2Thickness/2.0)),
                          niobiumBForAluminum2Logical,
                          "niobiumBForAluminum2Physical",worldLogical,false,0);

    //Now our final layer below aluminum2: aluminum3 (with daughter Nb)
    G4VSolid* aluminum3Solid
      = new G4Box("aluminum3Solid",0.1*cm,1.0*cm,dp_aluminum3Thickness/2.0);
    G4LogicalVolume* aluminum3Logical
      = new G4LogicalVolume(aluminum3Solid,fAluminum,"aluminum3Logical");
    G4VPhysicalVolume* aluminum3Physical
      = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,(-dp_siThickness/2.0 - dp_geThickness - dp_aluminum1Thickness - dp_aluminum2Thickness - dp_aluminum3Thickness/2.0)),
                          aluminum3Logical,"aluminum3Physical", worldLogical,
                          false,0);

    G4VSolid* niobiumForAluminum3Solid
      = new G4Box("niobiumForAluminum3Solid",0.1*cm,0.001*cm,
                  dp_aluminum3Thickness/2.0);
    G4LogicalVolume* niobiumForAluminum3Logical
      = new G4LogicalVolume(niobiumForAluminum3Solid,fNiobium,
                            "niobiumForAluminum3Logical");
    G4VPhysicalVolume* niobiumForAluminum3Physical
      = new G4PVPlacement(0,G4ThreeVector(),niobiumForAluminum3Logical,
                          "niobiumForAluminum3Physical",aluminum3Logical,false,
                          0);




    
    //Define lattices next
    G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();
    G4LatticeLogical* SiLogical = LM->LoadLattice(fSilicon, "Si");
    G4LatticeLogical* GeLogical = LM->LoadLattice(fGermanium, "Ge");
    G4LatticeLogical* AlLogical = LM->LoadLattice(fAluminum, "Al");
    G4LatticeLogical* NbLogical = LM->LoadLattice(fNiobium, "Nb");

    // G4LatticePhysical assigns G4LatticeLogical a physical orientation
    //For Ge, making Si for now because Ge has some reflection issues that
    //I'll wait on NT's push to fix
    G4LatticePhysical* SiPhysical = new G4LatticePhysical(SiLogical);
    G4LatticePhysical* GePhysical = new G4LatticePhysical(GeLogical);
    G4LatticePhysical* Al1Physical =
      new G4LatticePhysical(AlLogical,dp_polycryElScatMFP_Al,
                            dp_scDelta0_Al, dp_scTeff_Al,
                            dp_scDn_Al, dp_scTauQPTrap_Al);
    G4LatticePhysical* Al2Physical =
      new G4LatticePhysical(AlLogical,dp_polycryElScatMFP_Al,
                            dp_scDelta0_Al, dp_scTeff_Al,
                            dp_scDn_Al, dp_scTauQPTrap_Al);
    G4LatticePhysical* Al3Physical =
      new G4LatticePhysical(AlLogical,dp_polycryElScatMFP_Al,
                            dp_scDelta0_Al, dp_scTeff_Al,
                            dp_scDn_Al, dp_scTauQPTrap_Al);
    G4LatticePhysical* NbAForAl2Physical =
      new G4LatticePhysical(NbLogical,dp_polycryElScatMFP_Nb,
                            dp_scDelta0_Nb, dp_scTeff_Nb,
                            dp_scDn_Nb, dp_scTauQPTrap_Nb);
    G4LatticePhysical* NbBForAl2Physical =
      new G4LatticePhysical(NbLogical,dp_polycryElScatMFP_Nb,
                            dp_scDelta0_Nb, dp_scTeff_Nb,
                            dp_scDn_Nb, dp_scTauQPTrap_Nb);
    G4LatticePhysical* NbForAl3Physical =
      new G4LatticePhysical(NbLogical,dp_polycryElScatMFP_Nb,
                            dp_scDelta0_Nb, dp_scTeff_Nb,
                            dp_scDn_Nb, dp_scTauQPTrap_Nb);
    SiPhysical->SetMillerOrientation(1,0,0);
    GePhysical->SetMillerOrientation(1,0,0);
    Al1Physical->SetMillerOrientation(1,0,0);
    Al2Physical->SetMillerOrientation(1,0,0);
    Al3Physical->SetMillerOrientation(1,0,0);
    NbAForAl2Physical->SetMillerOrientation(1,0,0);
    NbBForAl2Physical->SetMillerOrientation(1,0,0);
    NbForAl3Physical->SetMillerOrientation(1,0,0);
    LM->RegisterLattice(siliconPhysical, SiPhysical);
    LM->RegisterLattice(germaniumPhysical, GePhysical);
    LM->RegisterLattice(aluminum1Physical, Al1Physical);
    LM->RegisterLattice(aluminum2Physical, Al2Physical);
    LM->RegisterLattice(aluminum3Physical, Al3Physical);
    LM->RegisterLattice(niobiumAForAluminum2Physical, NbAForAl2Physical);
    LM->RegisterLattice(niobiumBForAluminum2Physical, NbBForAl2Physical);
    LM->RegisterLattice(niobiumForAluminum3Physical, NbForAl3Physical);

    //Now do Border surface properties
    if (!fConstructed) {
      G4cout << "Constructing borders." << G4endl;
      
      const G4double GHz = 1e9 * hertz; 
      
      //the following coefficients and cutoff values are not well-motivated
      //the code below is used only to demonstrate how to set these values.
      const std::vector<G4double> anhCoeffs = {0, 0, 0, 0, 0, 1.51e-14};
      const std::vector<G4double> diffCoeffs = {};
      const std::vector<G4double> specCoeffs = {};
    
      const G4double anhCutoff = 520., reflCutoff = 350.;   // Units external

      //Boundary between any surface and vacuum: fully reflective
      fVacSurfProp = new G4CMPSurfaceProperty("VacSurf",
                                              0.0, 1.0, 0.0, 0.0,
                                              0.0, 1.0, 0.0, 0.0,
                                              0.0, 1.0);
      fVacSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
                                            diffCoeffs, specCoeffs, GHz, GHz,
                                            GHz);

      //Silicon-germanium boundary: 50% reflection for phonons, absorbant for
      //QPs (bc they aren't realizable here)
      fSiGeSurfProp = new G4CMPSurfaceProperty("SiGeSurf",
                                               0.0, 1.0, 0.0, 0.0,
                                               0.0, 0.5, 0.0, 0.0,
                                               1.0, 0.0);
      fSiGeSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
                                             diffCoeffs, specCoeffs, GHz, GHz,
                                             GHz);
      
      //Germanium-Aluminum1 boundary, 100% reflection, to isolate SC from non-SC
      fGeAl1SurfProp = new G4CMPSurfaceProperty("GeAl1Surf",
                                                0.0, 1.0, 0.0, 0.0,
                                                0.0, 1.0, 0.0, 0.0,
                                                0.0, 1.0);
      fGeAl1SurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
                                              diffCoeffs, specCoeffs, GHz, GHz,
                                              GHz);

      //Aluminum1-Aluminum2 boundary, from above: 100% absorption, to kill
      //phonons in quasielastic scattering
      fAl1Al2FromAboveSurfProp = new G4CMPSurfaceProperty("Al1Al2FromAboveSurf",
                                                          0.0, 1.0, 0.0, 0.0,
                                                          1.0, 0.0, 0.0, 0.0,
                                                          0.0, 1.0);
      fAl1Al2FromAboveSurfProp->AddScatteringProperties(anhCutoff, reflCutoff,
                                                        anhCoeffs,diffCoeffs,
                                                        specCoeffs, GHz, GHz,
                                                        GHz);
      
      //Aluminum1-Aluminum2 boundary, from below: 100% reflection, to isolate
      //Aluminum 2 region
      fAl1Al2FromBelowSurfProp = new G4CMPSurfaceProperty("Al1Al2FromBelowSurf",
                                                          0.0, 1.0, 0.0, 0.0,
                                                          0.0, 1.0, 0.0, 0.0,
                                                          0.0, 1.0);
      fAl1Al2FromBelowSurfProp->AddScatteringProperties(anhCutoff, reflCutoff,
                                                        anhCoeffs,diffCoeffs,
                                                        specCoeffs, GHz, GHz,
                                                        GHz);
      
      //Aluminum2-Niobium boundary: 100% QP absorption, 100% phonon reflection,
      //to do TOFP plots
      fAl2NbSurfProp = new G4CMPSurfaceProperty("Al2NbSurf",
                                                0.0, 1.0, 0.0, 0.0,
                                                0.0, 1.0, 0.0, 0.0,
                                                1.0, 0.0);
      fAl2NbSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
                                              diffCoeffs, specCoeffs, GHz, GHz,
                                              GHz);

      //Aluminum2-Aluminum3 boundary: 100% reflection, to isolate Aluminum 1
      //region
      fAl2Al3SurfProp = new G4CMPSurfaceProperty("Al2Al3Surf",
                                                 0.0, 1.0, 0.0, 0.0,
                                                 0.0, 1.0, 0.0, 0.0,
                                                 0.0, 1.0);
      fAl2Al3SurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
                                               diffCoeffs, specCoeffs, GHz, GHz,
                                               GHz);

      //Aluminum3-Niobium boundary: 100% QP and phonon transmission, to enable
      //proper gap studies
      fAl3NbSurfProp = new G4CMPSurfaceProperty("Al3NbSurf",
                                                0.0, 1.0, 0.0, 0.0,
                                                0.0, 0.0, 0.0, 0.0,
                                                0.0, 0.0);
      fAl3NbSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
                                              diffCoeffs, specCoeffs, GHz, GHz,
                                              GHz);

      
    }


    //Now we do the actual logical border surface creation
    //Silicon
    new G4CMPLogicalBorderSurface("SiVac",siliconPhysical,fWorldPhys,
                                  fVacSurfProp);
    new G4CMPLogicalBorderSurface("VacSi",fWorldPhys,siliconPhysical,
                                  fVacSurfProp);    
    new G4CMPLogicalBorderSurface("SiGe",siliconPhysical,germaniumPhysical,
                                  fSiGeSurfProp);
    new G4CMPLogicalBorderSurface("GeSi",germaniumPhysical,siliconPhysical,
                                  fSiGeSurfProp);

    //Germanium
    new G4CMPLogicalBorderSurface("GeVac",germaniumPhysical,fWorldPhys,
                                  fVacSurfProp);
    new G4CMPLogicalBorderSurface("VacGe",fWorldPhys,germaniumPhysical,
                                  fVacSurfProp);    
    new G4CMPLogicalBorderSurface("GeAl1",germaniumPhysical,aluminum1Physical,
                                  fGeAl1SurfProp);
    new G4CMPLogicalBorderSurface("Al1Ge",aluminum1Physical,germaniumPhysical,
                                  fGeAl1SurfProp);

    //Aluminum1
    new G4CMPLogicalBorderSurface("Al1Vac",aluminum1Physical,fWorldPhys,
                                  fVacSurfProp);
    new G4CMPLogicalBorderSurface("VacAl1",fWorldPhys,aluminum1Physical,
                                  fVacSurfProp);
    new G4CMPLogicalBorderSurface("Al1Al2FromAbove",aluminum1Physical,
                                  aluminum2Physical,fAl1Al2FromAboveSurfProp);
    new G4CMPLogicalBorderSurface("Al2Al1FromBelow",aluminum2Physical,
                                  aluminum1Physical,fAl1Al2FromBelowSurfProp);

    //Aluminum2
    new G4CMPLogicalBorderSurface("Al2Vac",aluminum2Physical,fWorldPhys,
                                  fVacSurfProp);
    new G4CMPLogicalBorderSurface("VacAl2",fWorldPhys,aluminum2Physical,
                                  fVacSurfProp);
    new G4CMPLogicalBorderSurface("Al2NbA",aluminum2Physical,
                                  niobiumAForAluminum2Physical,fAl2NbSurfProp);
    new G4CMPLogicalBorderSurface("NbAAl2",niobiumAForAluminum2Physical,
                                  aluminum2Physical,fAl2NbSurfProp);
    new G4CMPLogicalBorderSurface("Al2NbB",aluminum2Physical,
                                  niobiumBForAluminum2Physical,fAl2NbSurfProp);
    new G4CMPLogicalBorderSurface("NbBAl2",niobiumBForAluminum2Physical,
                                  aluminum2Physical,fAl2NbSurfProp);
    new G4CMPLogicalBorderSurface("Al2Al3",aluminum2Physical,aluminum3Physical,
                                  fAl2Al3SurfProp);
    new G4CMPLogicalBorderSurface("Al3Al2",aluminum3Physical,aluminum2Physical,
                                  fAl2Al3SurfProp);    
    new G4CMPLogicalBorderSurface("Al2NbLower",aluminum2Physical,
                                  niobiumForAluminum3Physical,fAl2Al3SurfProp);
    new G4CMPLogicalBorderSurface("NbLowerAl2",niobiumForAluminum3Physical,
                                  aluminum2Physical,fAl2Al3SurfProp);

    
    //Aluminum2's niobium ends
    new G4CMPLogicalBorderSurface("NbAVac",niobiumAForAluminum2Physical,
                                  fWorldPhys,fVacSurfProp);
    new G4CMPLogicalBorderSurface("VacNbA",fWorldPhys,
                                  niobiumAForAluminum2Physical,fVacSurfProp);
    new G4CMPLogicalBorderSurface("NbBVac",niobiumBForAluminum2Physical,
                                  fWorldPhys,fVacSurfProp);
    new G4CMPLogicalBorderSurface("VacNbB",fWorldPhys,
                                  niobiumBForAluminum2Physical,fVacSurfProp);
    
    //Aluminum3
    new G4CMPLogicalBorderSurface("Al3Vac",aluminum3Physical,fWorldPhys,
                                  fVacSurfProp);
    new G4CMPLogicalBorderSurface("VacAl3",fWorldPhys,aluminum3Physical,
                                  fVacSurfProp);
    new G4CMPLogicalBorderSurface("NbAl3",niobiumForAluminum3Physical,
                                  aluminum3Physical,fAl3NbSurfProp);
    new G4CMPLogicalBorderSurface("Al3Nb",aluminum3Physical,
                                  niobiumForAluminum3Physical,fAl3NbSurfProp);

    //Aluminum3's niobium center
    new G4CMPLogicalBorderSurface("NbVac",niobiumForAluminum3Physical,
                                  fWorldPhys,fVacSurfProp);
    new G4CMPLogicalBorderSurface("VacNb",fWorldPhys,
                                  niobiumForAluminum3Physical,fVacSurfProp);




    
    //Set visualization attributes
    worldLogical->SetVisAttributes(G4VisAttributes::Invisible);
    G4VisAttributes* vacuumBoxVisAtt
      = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    
    vacuumBoxVisAtt->SetVisibility(true);    
    G4VisAttributes* siliconBoxVisAtt
      = new G4VisAttributes(G4Colour(0.5,0.5,0.5,0.5));
    
    siliconLogical->SetVisAttributes(siliconBoxVisAtt);
    siliconBoxVisAtt->SetVisibility(true);    
    G4VisAttributes* germaniumBoxVisAtt
      = new G4VisAttributes(G4Colour(0.5,0.5,0.5,0.5));
    
    germaniumLogical->SetVisAttributes(germaniumBoxVisAtt);
    germaniumBoxVisAtt->SetVisibility(true);    
    G4VisAttributes* aluminumBoxVisAtt
      = new G4VisAttributes(G4Colour(1.0,0.0,1.0,0.5));
    
    aluminum1Logical->SetVisAttributes(aluminumBoxVisAtt);
    aluminum2Logical->SetVisAttributes(aluminumBoxVisAtt);
    aluminum3Logical->SetVisAttributes(aluminumBoxVisAtt);
    aluminumBoxVisAtt->SetVisibility(true);    
    G4VisAttributes* niobiumBoxVisAtt
      = new G4VisAttributes(G4Colour(0.0,1.0,1.0,0.5));
    
    niobiumAForAluminum2Logical->SetVisAttributes(niobiumBoxVisAtt);
    niobiumBForAluminum2Logical->SetVisAttributes(niobiumBoxVisAtt);
    niobiumForAluminum3Logical->SetVisAttributes(niobiumBoxVisAtt);
    niobiumBoxVisAtt->SetVisibility(true);
    
    
  } else if (geometryID == 2) {
    
    //Start by defining interface properties (since these are needed by
    //classes we instantiate objects of here)
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
                                              diffCoeffs, specCoeffs, GHz, GHz,
                                              GHz);
      fBorderContainer.emplace("SiAl",fSiAlInterface);

      //For the the interface of the Si and the world
      fSiVacInterface = new G4CMPSurfaceProperty("SiVacSurf",
                                                 0.0, 1.0, 0.0, 0.0,
                                                 0.0, 1.0, 0.0, 0.0,
                                                 0.0, 1.0);
      fSiVacInterface->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
                                               diffCoeffs, specCoeffs, GHz, GHz,
                                               GHz);
      fBorderContainer.emplace("SiVac",fSiVacInterface);

      //For the the interface of the Si and the Cu
      fSiCuInterface = new G4CMPSurfaceProperty("SiCuSurf",
                                                0.0, 1.0, 0.0, 0.0,
                                                0.0, 0.0, 0.0, 0.0,
                                                0.0, 1.0);
      fSiCuInterface->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
                                              diffCoeffs, specCoeffs, GHz, GHz,
                                              GHz);
      fBorderContainer.emplace("SiCu",fSiCuInterface);      

      //For the the interface of the Cu and the Vac
      fCuVacInterface = new G4CMPSurfaceProperty("CuVacSurf",
                                                 0.0, 1.0, 0.0, 0.0,
                                                 1.0, 1.0, 0.0, 0.0,
                                                 0.0, 1.0);
      fCuVacInterface->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
                                               diffCoeffs, specCoeffs, GHz, GHz,
                                               GHz);
      fBorderContainer.emplace("CuVac",fCuVacInterface);      

      
      //For the the interface of the Al and world
      fAlVacInterface = new G4CMPSurfaceProperty("AlVacSurf",
                                                 0.0, 1.0, 0.0, 0.0,
                                                 0.0, 1.0, 0.0, 0.0,
                                                 0.0, 1.0);
      fAlVacInterface->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
                                               diffCoeffs, specCoeffs, GHz, GHz,
                                               GHz);
      fBorderContainer.emplace("AlVac",fAlVacInterface);      

      //For the the interface of the Al and Al
      fAlAlInterface = new G4CMPSurfaceProperty("AlAlSurf",
                                                0.0, 1.0, 0.0, 0.0,
                                                0.0, 0.0, 0.0, 0.0,
                                                0.0, 0.0);
      fAlAlInterface->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
                                              diffCoeffs, specCoeffs, GHz, GHz,
                                              GHz);
      fBorderContainer.emplace("AlAl",fAlAlInterface);

      //For the the interface of the Vac and Vac
      fVacVacInterface = new G4CMPSurfaceProperty("VacVacSurf",
                                                  0.0, 1.0, 0.0, 0.0,
                                                  0.0, 1.0, 0.0, 0.0,
                                                  0.0, 1.0);
      fVacVacInterface->AddScatteringProperties(anhCutoff, reflCutoff,
                                                anhCoeffs,diffCoeffs,
                                                specCoeffs,GHz,GHz,GHz);
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
       
    // Now we start constructing the various components and their interfaces  
    bool checkOverlaps = true;

    //First, set up the qubit chip substrate. 
    G4Box * solid_siliconChip = new G4Box("QubitChip_solid",
                                          0.5*dp_siliconChipDimX,
                                          0.5*dp_siliconChipDimY,
                                          0.5*dp_siliconChipDimZ);
  
    //Now attribute a physical material to the chip
    G4LogicalVolume * log_siliconChip = new G4LogicalVolume(solid_siliconChip,
                                                            fSilicon,
                                                            "SiliconChip_log");
    
    //Now, create a physical volume and G4PVPlacement for storing as the
    //final output
    G4ThreeVector siliconChipTranslate(0,0,0.5*(dp_housingDimZ - dp_siliconChipDimZ) + dp_eps); 
    G4VPhysicalVolume * phys_siliconChip
      = new G4PVPlacement(0,siliconChipTranslate,log_siliconChip,"SiliconChip",
                          worldLogical,false,0,checkOverlaps);
    
    G4VisAttributes* siliconChipVisAtt
      = new G4VisAttributes(G4Colour(0.5,0.5,0.5));
    siliconChipVisAtt->SetVisibility(true);
    log_siliconChip->SetVisAttributes(siliconChipVisAtt);



    // G4LatticePhysical assigns G4LatticeLogical a physical orientation
    G4LatticePhysical* phys_siliconLattice
      = new G4LatticePhysical(log_siliconLattice);
    phys_siliconLattice->SetMillerOrientation(1,0,0); 
    LM->RegisterLattice(phys_siliconChip,phys_siliconLattice);

    //Set up border surfaces
    G4CMPLogicalBorderSurface * border_siliconChip_world
      = new G4CMPLogicalBorderSurface("border_siliconChip_world",
                                      phys_siliconChip, fWorldPhys,
                                      fSiVacInterface);

    
    //If desired, set up the copper qubit housing
    if (dp_useQubitHousing) {           
      ValidationQubitHousing * qubitHousing
        = new ValidationQubitHousing(0,G4ThreeVector(0,0,0),"QubitHousing",
                                     worldLogical,false,0,checkOverlaps);
      G4LogicalVolume * log_qubitHousing = qubitHousing->GetLogicalVolume();
      G4VPhysicalVolume * phys_qubitHousing = qubitHousing->GetPhysicalVolume();

      //Set up lattice properties for copper housing
      G4LatticePhysical* phys_copperLattice
        = new G4LatticePhysical(log_copperLattice);
      phys_copperLattice->SetMillerOrientation(1,0,0); 
      LM->RegisterLattice(phys_qubitHousing,phys_copperLattice);


      
      //Set up the logical border surface
      G4CMPLogicalBorderSurface * border_siliconChip_qubitHousing
        = new G4CMPLogicalBorderSurface("border_siliconChip_qubitHousing",
                                        phys_siliconChip, phys_qubitHousing,
                                        fSiCuInterface);
      G4CMPLogicalBorderSurface * border_qubitHousing_siliconChip
        = new G4CMPLogicalBorderSurface("border_qubitHousing_siliconChip",
                                        phys_qubitHousing, phys_siliconChip,
                                        fSiCuInterface);
      G4CMPLogicalBorderSurface * border_qubitHousing_world
        = new G4CMPLogicalBorderSurface("border_qubitHousing_world",
                                        phys_qubitHousing, fWorldPhys,
                                        fCuVacInterface);

    }

    
    //Now set up the ground plane, in which the transmission line, resonators, and qubits will be located.
    if (dp_useGroundPlane) {    
    
      G4Box * solid_groundPlane
        = new G4Box("GroundPlane_solid",0.5*dp_groundPlaneDimX,
                    0.5*dp_groundPlaneDimY,0.5*dp_groundPlaneDimZ);
    
    
      //Now attribute a physical material to the chip
      G4LogicalVolume * log_groundPlane
        = new G4LogicalVolume(solid_groundPlane,fAluminum,"GroundPlane_log");
    
    
      //Now, create a physical volume and G4PVPlacement for storing as the
      //final output
      G4ThreeVector groundPlaneTranslate(0,0,0.5*(dp_housingDimZ) + dp_eps + dp_groundPlaneDimZ*0.5);
      G4VPhysicalVolume * phys_groundPlane
        = new G4PVPlacement(0,groundPlaneTranslate,log_groundPlane,
                            "GroundPlane",worldLogical,false,0,checkOverlaps);
      
      G4VisAttributes* groundPlaneVisAtt
        = new G4VisAttributes(G4Colour(0.0,1.0,1.0,0.5));
      groundPlaneVisAtt->SetVisibility(true);
      log_groundPlane->SetVisAttributes(groundPlaneVisAtt);

      G4LatticePhysical* phys_groundPlaneLattice
        = new G4LatticePhysical(log_aluminumLattice,dp_polycryElScatMFP_Al,
                                dp_scDelta0_Al, dp_scTeff_Al,
                                dp_scDn_Al, dp_scTauQPTrap_Al);
      phys_groundPlaneLattice->SetMillerOrientation(1,0,0);
      LM->RegisterLattice(phys_groundPlane,phys_groundPlaneLattice);

      
      
      //Set up the logical border surface
      G4CMPLogicalBorderSurface * border_siliconChip_groundPlane
        = new G4CMPLogicalBorderSurface("border_siliconChip_groundPlane",
                                        phys_siliconChip, phys_groundPlane,
                                        fSiAlInterface);
      G4CMPLogicalBorderSurface * border_groundPlane_siliconChip
        = new G4CMPLogicalBorderSurface("border_siliconChip_groundPlane",
                                        phys_groundPlane, phys_siliconChip,
                                        fSiAlInterface);
      G4CMPLogicalBorderSurface * border_world_groundPlane
        = new G4CMPLogicalBorderSurface("border_world_groundPlane", fWorldPhys,
                                        phys_groundPlane, fAlVacInterface);
      G4CMPLogicalBorderSurface * border_groundPlane_world
        = new G4CMPLogicalBorderSurface("border_groundPlane_world",
                                        phys_groundPlane, fWorldPhys,
                                        fAlVacInterface);


    
    
      //Now set up the transmission line
      if (dp_useTransmissionLine) {
	
        //Since it's within the ground plane exactly;
        //0.5*(dp_housingDimZ) + dp_eps + dp_groundPlaneDimZ*0.5 ); 
        G4ThreeVector transmissionLineTranslate(0,0,0.0);
        ValidationTransmissionLine * tLine
          = new ValidationTransmissionLine(0,transmissionLineTranslate,
                                           "TransmissionLine",log_groundPlane,
                                           false,0,LM,fLogicalLatticeContainer,
                                           fBorderContainer,checkOverlaps);
        G4LogicalVolume * log_tLine = tLine->GetLogicalVolume();
        G4VPhysicalVolume * phys_tLine = tLine->GetPhysicalVolume();



        //Now, if we're using the chip and ground plane AND the transmission
        //line. This gets a bit hairy, since the transmission line is
        //composite of both Nb and vacuum. So we'll access the list of
        //physical objects present in it and link those one-by-one to the
        //silicon chip and to the world
        for (int iSubVol = 0; iSubVol < tLine->GetListOfAllFundamentalSubVolumes().size(); ++iSubVol) {

          std::tuple<std::string,G4String,G4VPhysicalVolume*> theTLTuple
            = tLine->GetListOfAllFundamentalSubVolumes()[iSubVol];

          std::cout << "TLine sub volume names (to be used for boundaries): "
                    << std::get<1>(theTLTuple)
                    << " with material "
                    << std::get<0>(theTLTuple)
                    << std::endl;

          //Set the si/vac and si/Al interfaces to be symmetric in both
          //dimensions
          std::string tempName1a = "border_siliconChip_"
            + std::get<1>(theTLTuple);
          std::string tempName1b = "border_" + std::get<1>(theTLTuple)
            + "_siliconChip";
	  
          if (std::get<0>(theTLTuple).find("Vacuum") != std::string::npos) {

            G4CMPLogicalBorderSurface *
              border1_siliconChip_transmissionLineEmpty
              = new G4CMPLogicalBorderSurface(tempName1a, phys_siliconChip,
                                              std::get<2>(theTLTuple),
                                              fSiVacInterface);
	    
            G4CMPLogicalBorderSurface *
              border2_siliconChip_transmissionLineEmpty
              = new G4CMPLogicalBorderSurface(tempName1b,
                                              std::get<2>(theTLTuple),
                                              phys_siliconChip,
                                              fSiVacInterface);
          }
	  
          if (std::get<0>(theTLTuple).find("Aluminum") != std::string::npos) {
	    
            G4CMPLogicalBorderSurface *
              border1_siliconChip_transmissionLineConductor
              = new G4CMPLogicalBorderSurface(tempName1a, phys_siliconChip,
                                              std::get<2>(theTLTuple),
                                              fSiAlInterface);

            G4CMPLogicalBorderSurface *
              border2_siliconChip_transmissionLineConductor
              = new G4CMPLogicalBorderSurface(tempName1b,
                                              std::get<2>(theTLTuple),
                                              phys_siliconChip,
                                              fSiAlInterface);
          }


          //Set the world/vac and world/Al interfaces to be symmetric in both
          //dimensions
          std::string tempName2a = "border_world_" + std::get<1>(theTLTuple);
          std::string tempName2b = "border_" + std::get<1>(theTLTuple) +
            "_world";
	  
          if (std::get<0>(theTLTuple).find("Vacuum") != std::string::npos) {

            G4CMPLogicalBorderSurface *
              border1_world_transmissionLineEmpty
              = new G4CMPLogicalBorderSurface(tempName2a, fWorldPhys,
                                              std::get<2>(theTLTuple),
                                              fVacVacInterface);
            G4CMPLogicalBorderSurface *
              border2_world_transmissionLineEmpty
              = new G4CMPLogicalBorderSurface(tempName2b,
                                              std::get<2>(theTLTuple),
                                              fWorldPhys,fVacVacInterface);
          }
          if (std::get<0>(theTLTuple).find("Aluminum") != std::string::npos) {
	    
            G4CMPLogicalBorderSurface *
              border1_world_transmissionLineConductor
              = new G4CMPLogicalBorderSurface(tempName2a, fWorldPhys,
                                              std::get<2>(theTLTuple),
                                              fAlVacInterface);
            G4CMPLogicalBorderSurface *
              border2_world_transmissionLineConductor
              = new G4CMPLogicalBorderSurface(tempName2b,
                                              std::get<2>(theTLTuple),
                                              fWorldPhys,fAlVacInterface);
          }
	  
	  
          //Space here for linking to GROUND PLANE. Only need to link vacuum
          //ones since they're fully enclosing TL
          std::string tempName3a = "border_groundPlane_" +
            std::get<1>(theTLTuple);
          std::string tempName3b = "border_" + std::get<1>(theTLTuple) +
            "_groundPlane";	  
          if (std::get<0>(theTLTuple).find("Vacuum") != std::string::npos) {
            G4CMPLogicalBorderSurface *
              border1_groundPlane_transmissionLineEmpty
              = new G4CMPLogicalBorderSurface(tempName3a, phys_groundPlane,
                                              std::get<2>(theTLTuple),
                                              fAlVacInterface);
	    
            G4CMPLogicalBorderSurface *
              border2_groundPlane_transmissionLineEmpty
              = new G4CMPLogicalBorderSurface(tempName3b,
                                              std::get<2>(theTLTuple),
                                              phys_groundPlane,
                                              fAlVacInterface);
          }
        }
      }
    
  
      //Now set up a set of 6 resonator assemblies
      if (dp_useResonatorAssembly) {
        int nR = 6;
        for (int iR = 0; iR < nR; ++iR) {
      
          //First, get the translation vector for the resonator assembly
          //For the top three, don't do a rotation. For the bottom three, do
          G4ThreeVector resonatorAssemblyTranslate(0,0,0);
          G4RotationMatrix * rotAssembly = 0;
          if (iR <= 2) {
            resonatorAssemblyTranslate
              = G4ThreeVector(dp_resonatorLateralSpacing*(iR-1)
                              +dp_centralResonatorOffsetX,
                              0.5 * dp_resonatorAssemblyBaseAlDimY
                              + 0.5 * dp_transmissionLineCavityFullWidth,
                              0.0);
            rotAssembly = 0;
          } else {

            //Negative offset because qubit is mirrored on underside
            resonatorAssemblyTranslate
              = G4ThreeVector(dp_resonatorLateralSpacing*(iR-4)-
                              dp_centralResonatorOffsetX, 
                              -1*(0.5 * dp_resonatorAssemblyBaseAlDimY +
                                  0.5 * dp_transmissionLineCavityFullWidth),
                              0.0);
            rotAssembly = new G4RotationMatrix();
            rotAssembly->rotateZ(180*deg);
          }
	
          char name[400];
          sprintf(name,"ResonatorAssembly_%d",iR);
          G4String resonatorAssemblyName(name);
          ValidationResonatorAssembly * resonatorAssembly
            = new ValidationResonatorAssembly(rotAssembly,
                                              resonatorAssemblyTranslate,
                                              resonatorAssemblyName,
                                              log_groundPlane,false,
                                              0,LM,fLogicalLatticeContainer,
                                              fBorderContainer,checkOverlaps);

          G4LogicalVolume * log_resonatorAssembly
            = resonatorAssembly->GetLogicalVolume();
          G4VPhysicalVolume * phys_resonatorAssembly
            = resonatorAssembly->GetPhysicalVolume();

          //Do the logical border creation now
          for (int iSubVol = 0; iSubVol < resonatorAssembly->GetListOfAllFundamentalSubVolumes().size(); ++iSubVol) {

            std::tuple<std::string,G4String,G4VPhysicalVolume*> theResTuple
              = resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol];
	    
            std::cout << "Resonator assuembly sub volume names (to be used for "
                      << "boundaries): "
                      << std::get<1>(theResTuple) << " with material "
                      << std::get<0>(theResTuple) << std::endl;

	    
            //Set the chip/vacuum interfaces
            if (std::get<0>(theResTuple).find("Vacuum") != std::string::npos) {
              std::string tempName1 = "border_siliconChip_"
                + std::get<1>(theResTuple);
              std::string tempName2 = "border_"
                + std::get<1>(theResTuple) + "_siliconChip";
              G4CMPLogicalBorderSurface *
                border_siliconChip_resonatorAssemblyEmpty
                = new G4CMPLogicalBorderSurface(tempName1,
                                                std::get<2>(theResTuple),
                                                phys_siliconChip,
                                                fSiVacInterface);
              G4cout << "Setting interface for " << tempName1 << " to fSiVacInterface" << G4endl;
	      
              G4CMPLogicalBorderSurface *
                border_resonatorAssemblyEmpty_siliconChip
                = new G4CMPLogicalBorderSurface(tempName2,
                                                phys_siliconChip,
                                                std::get<2>(theResTuple),
                                                fSiVacInterface);
              G4cout << "Setting interface for " << tempName2 << " to fSiVacInterface" << G4endl;
            }

            //Set the chip/aluminum interfaces
            if (std::get<0>(theResTuple).find("Aluminum") !=std::string::npos) {
              std::string tempName1 = "border_siliconChip_" +
                std::get<1>(theResTuple);
              std::string tempName2 = "border_" + std::get<1>(theResTuple) +
                "_siliconChip";
              G4CMPLogicalBorderSurface *
                border_siliconChip_resonatorAssemblyConductor
                = new G4CMPLogicalBorderSurface(tempName1, phys_siliconChip,
                                                std::get<2>(theResTuple),
                                                fSiAlInterface);
              G4cout << "Setting interface for " << tempName1 << " to fSiAlInterface" << G4endl;
	      
              G4CMPLogicalBorderSurface *
                border_resonatorAssemblyConductor_siliconChip
                = new G4CMPLogicalBorderSurface(tempName2,
                                                std::get<2>(theResTuple),
                                                phys_siliconChip,
                                                fSiAlInterface);
              G4cout << "Setting interface for " << tempName2 << " to fSiAlInterface" << G4endl;
            }
	    
            //Set the world/vacuum interfaces (probably not necessary but whatever, better to have everything well-defined
            if (std::get<0>(theResTuple).find("Vacuum") != std::string::npos) {
              std::string tempName1 = "border_world_" +
                std::get<1>(theResTuple);
              std::string tempName2 = "border_" + std::get<1>(theResTuple) +
                "_world";
              G4CMPLogicalBorderSurface *
                border_world_resonatorAssemblyEmpty
                = new G4CMPLogicalBorderSurface(tempName1, fWorldPhys,
                                                std::get<2>(theResTuple),
                                                fVacVacInterface);
              G4cout << "Setting interface for " << tempName1 << " to fVacVacInterface" << G4endl;
	      
              G4CMPLogicalBorderSurface *
                border_resonatorAssemblyEmpty_world
                = new G4CMPLogicalBorderSurface(tempName2,
                                                std::get<2>(theResTuple),
                                                fWorldPhys, fVacVacInterface);
              G4cout << "Setting interface for " << tempName2 << " to fVacVacInterface" << G4endl;
            }

            //Set the world/aluminum interfaces
            if (std::get<0>(theResTuple).find("Aluminum") !=std::string::npos) {
              std::string tempName1 = "border_world_" +
                std::get<1>(theResTuple);
              std::string tempName2 = "border_" + std::get<1>(theResTuple) +
                "_world";
              G4CMPLogicalBorderSurface *
                border_world_resonatorAssemblyConductor
                = new G4CMPLogicalBorderSurface(tempName1, fWorldPhys,
                                                std::get<2>(theResTuple),
                                                fAlVacInterface);
              G4cout << "Setting interface for " << tempName1 << " to fAlVacInterface" << G4endl;
	      
              G4CMPLogicalBorderSurface *
                border_resonatorAssemblyConductor_world
                = new G4CMPLogicalBorderSurface(tempName2,
                                                std::get<2>(theResTuple),
                                                fWorldPhys, fAlVacInterface);
              G4cout << "Setting interface for " << tempName2 << " to fAlVacInterface" << G4endl;
            }

            //Set the TLcouplerConductor interface with the ground plane
            if (std::get<1>(theResTuple).find("tlCouplingConductor") != std::string::npos) {
              std::string tempName1 = "border_groundPlane_" +
                std::get<1>(theResTuple);
              std::string tempName2 = "border_" + std::get<1>(theResTuple) +
                "_groundPlane";
              G4CMPLogicalBorderSurface *
                border_groundPlane_tlCouplingConductor
                = new G4CMPLogicalBorderSurface(tempName1, phys_groundPlane,
                                                std::get<2>(theResTuple),
                                                fAlAlInterface);
              G4cout << "Setting interface for " << tempName1 << " to fAlAlInterface" << G4endl;
              G4CMPLogicalBorderSurface *
                border_tlCouplingConductor_groundPlane
                = new G4CMPLogicalBorderSurface(tempName2,
                                                std::get<2>(theResTuple),
                                                phys_groundPlane,
                                                fAlAlInterface);
              G4cout << "Setting interface for " << tempName2 << " to fAlAlInterface" << G4endl;
            }
	    
	    
            //Set the TLCouplerEmpty interface with the ground plane
            if (std::get<1>(theResTuple).find("tlCouplingEmpty") != std::string::npos) {
              std::string tempName1 = "border_groundPlane_" +
                std::get<1>(theResTuple);
              std::string tempName2 = "border_" + std::get<1>(theResTuple) +
                "_groundPlane";
              G4CMPLogicalBorderSurface *
                border_groundPlane_tlCouplingEmpty
                = new G4CMPLogicalBorderSurface(tempName1, phys_groundPlane,
                                                std::get<2>(theResTuple),
                                                fAlVacInterface);
              G4cout << "Setting interface for " << tempName1 << " to fAlVacInterface" << G4endl;
	      
              G4CMPLogicalBorderSurface *
                border_tlCouplingEmpty_groundPlane
                = new G4CMPLogicalBorderSurface(tempName2,
                                                std::get<2>(theResTuple),
                                                phys_groundPlane,
                                                fAlVacInterface);
              G4cout << "Setting interface for " << tempName2 << " to fAlVacInterface" << G4endl;
	      
            }

            //Set the baseLayer/groundplane interface
            if (std::get<1>(theResTuple) == resonatorAssemblyName) {
              std::string tempName1 = "border_groundPlane_" +
                std::get<1>(theResTuple);
              std::string tempName2 = "border_" + std::get<1>(theResTuple) +
                "_groundPlane";
              G4CMPLogicalBorderSurface * border_groundPlane_baselayer
                = new G4CMPLogicalBorderSurface(tempName1, phys_groundPlane,
                                                std::get<2>(theResTuple),
                                                fAlAlInterface);
              G4cout << "---Test Setting interface for " << tempName1 << " to fAlAlInterface" << G4endl;
	      
              G4CMPLogicalBorderSurface * border_baselayer_groundPlane
                = new G4CMPLogicalBorderSurface(tempName2,
                                                std::get<2>(theResTuple),
                                                phys_groundPlane,
                                                fAlAlInterface);
              G4cout << "---Test Setting interface for " << tempName2 << " to fAlAlInterface" << G4endl;
            }
          }
        }
      }
    }
  } else {
    
  }
}





// Attach material properties and electrode/sensor handler to surface

void ValidationDetectorConstruction::
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


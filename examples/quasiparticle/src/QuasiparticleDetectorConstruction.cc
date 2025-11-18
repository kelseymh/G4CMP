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
#include "QuasiparticleResonatorAssembly.hh"
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


QuasiparticleDetectorConstruction::QuasiparticleDetectorConstruction()
  : fLiquidHelium(0), fSilicon(0), fAluminum(0), fTungsten(0), fCopper(0),
    fNiobium(0),
    fWorldPhys(0), topSurfProp(0), vacSurfProp(0), wallSurfProp(0),
    copTopSurfProp(0), copTopSurfProp2(0), topSurfProp2(0), alNbSurfProp(0),
    electrodeSensitivity(0), fConstructed(false) {;}


QuasiparticleDetectorConstruction::~QuasiparticleDetectorConstruction() {
  delete topSurfProp;
  delete topSurfProp2;
  delete vacSurfProp;
  delete wallSurfProp;
  delete copTopSurfProp;
  delete copTopSurfProp2;
  delete alNbSurfProp;
}

G4VPhysicalVolume* QuasiparticleDetectorConstruction::Construct() {
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

void QuasiparticleDetectorConstruction::DefineMaterials() { 
  G4NistManager* nistManager = G4NistManager::Instance();

  //"Liquid helium" is a historical relic that should be replaced soon.
  //Here we're using it as a standin for vacuum
  fLiquidHelium = nistManager->FindOrBuildMaterial("G4_AIR"); 
  fSilicon = nistManager->FindOrBuildMaterial("G4_Si");
  fAluminum = nistManager->FindOrBuildMaterial("G4_Al");
  fCopper = nistManager->FindOrBuildMaterial("G4_Cu");
  fTungsten = nistManager->FindOrBuildMaterial("G4_W");
  fNiobium = nistManager->FindOrBuildMaterial("G4_Nb");
}

void QuasiparticleDetectorConstruction::SetupGeometry() {
  //     
  // World
  //
  G4VSolid* worldSolid = new G4Box("World",16.*cm,16.*cm,16.*cm);
  G4LogicalVolume* worldLogical =
    new G4LogicalVolume(worldSolid,fLiquidHelium,"World");
  worldLogical->SetUserLimits(new G4UserLimits(100*mm, DBL_MAX, DBL_MAX, 0, 0));
  fWorldPhys = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"World",0,
                                 false,0);

  
  //Start by defining interface properties (since these are needed by classes
  //we instantiate objects of here)
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
                                              0.0, 1.0, 0.0, 0.0,
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
    fVacVacInterface->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
                                              diffCoeffs, specCoeffs, GHz, GHz,
                                              GHz);
    fBorderContainer.emplace("VacVac",fVacVacInterface);
      
  }

  //Also need to define logical lattices *here* now, since the logical lattice
  //container needs to be passed into some classes
  // G4LatticeManager gives physics processes access to lattices by volume
  G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();
  G4LatticeLogical* log_siliconLattice = LM->LoadLattice(fSilicon, "Si");
  G4LatticeLogical* log_aluminumLattice = LM->LoadLattice(fAluminum, "Al");
  G4LatticeLogical* log_copperLattice = LM->LoadLattice(fCopper, "Cu");
  fLogicalLatticeContainer.emplace("Silicon",log_siliconLattice);
  fLogicalLatticeContainer.emplace("Aluminum",log_aluminumLattice);
  fLogicalLatticeContainer.emplace("Copper",log_copperLattice);

    

    
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  // Now we start constructing the various components and their interfaces  
  bool checkOverlaps = true;


  //--------------------------------------------------------------------------
  //First, set up the qubit chip substrate. 
  G4Box * solid_siliconChip = new G4Box("QubitChip_solid",
                                        0.5*dp_siliconChipDimX,
                                        0.5*dp_siliconChipDimY,
                                        0.5*dp_siliconChipDimZ);
  
  //Now attribute a physical material to the chip
  G4LogicalVolume * log_siliconChip = new G4LogicalVolume(solid_siliconChip,
                                                          fSilicon,
                                                          "SiliconChip_log");
    
  //Now, create a physical volume and G4PVPlacement for storing as the final
  //output
  G4ThreeVector siliconChipTranslate(0,0,
                                     0.5*(dp_housingDimZ - dp_siliconChipDimZ)
                                     + dp_eps); 
  G4VPhysicalVolume * phys_siliconChip =
    new G4PVPlacement(0,siliconChipTranslate,log_siliconChip,"SiliconChip",
                      worldLogical,false,0,checkOverlaps);
    
  G4VisAttributes* siliconChipVisAtt =
    new G4VisAttributes(G4Colour(0.5,0.5,0.5));

  siliconChipVisAtt->SetVisibility(true);
  log_siliconChip->SetVisAttributes(siliconChipVisAtt);

  // G4LatticePhysical assigns G4LatticeLogical a physical orientation
  G4LatticePhysical* phys_siliconLattice =
    new G4LatticePhysical(log_siliconLattice);
  phys_siliconLattice->SetMillerOrientation(1,0,0); 
  LM->RegisterLattice(phys_siliconChip,phys_siliconLattice);

  //Set up border surfaces
  G4CMPLogicalBorderSurface * border_siliconChip_world =
    new G4CMPLogicalBorderSurface("border_siliconChip_world",phys_siliconChip,
                                  fWorldPhys, fSiVacInterface);

    



  //--------------------------------------------------------------------------
  //If desired, set up the copper qubit housing
  if (dp_useQubitHousing) {

    QuasiparticleQubitHousing * qubitHousing =
      new QuasiparticleQubitHousing(0,G4ThreeVector(0,0,0),"QubitHousing",
                                    worldLogical,false,0,checkOverlaps);
    G4LogicalVolume * log_qubitHousing = qubitHousing->GetLogicalVolume();
    G4VPhysicalVolume * phys_qubitHousing = qubitHousing->GetPhysicalVolume();

    //Set up lattice properties for copper housing. The miller orientation is
    //a little silly for deliberately polycrystalline materials, but let's do
    //it for kicks anyway
    G4LatticePhysical* phys_copperLattice =
      new G4LatticePhysical(log_copperLattice);
    phys_copperLattice->SetMillerOrientation(1,0,0);
    LM->RegisterLattice(phys_qubitHousing,phys_copperLattice);
      
    //Set up the logical border surfaces
    G4CMPLogicalBorderSurface * border_siliconChip_qubitHousing =
      new G4CMPLogicalBorderSurface("border_siliconChip_qubitHousing",
                                    phys_siliconChip, phys_qubitHousing,
                                    fSiCuInterface);
    G4CMPLogicalBorderSurface * border_qubitHousing_siliconChip =
      new G4CMPLogicalBorderSurface("border_qubitHousing_siliconChip",
                                    phys_qubitHousing, phys_siliconChip,
                                    fSiCuInterface);
    G4CMPLogicalBorderSurface * border_qubitHousing_world =
      new G4CMPLogicalBorderSurface("border_qubitHousing_world",
                                    phys_qubitHousing, fWorldPhys,
                                    fCuVacInterface);
  }
    
  //-----------------------------------------------------------------
  //Now set up the ground plane, in which the transmission line, resonators,
  //and qubits will be located.
  if (dp_useGroundPlane) {

    G4Box * solid_groundPlane =
      new G4Box("GroundPlane_solid",0.5*dp_groundPlaneDimX,
                0.5*dp_groundPlaneDimY,0.5*dp_groundPlaneDimZ);
    
    //Now attribute a physical material to the chip
    G4LogicalVolume * log_groundPlane =
      new G4LogicalVolume(solid_groundPlane,fNiobium,"GroundPlane_log");
        
    //Now, create a physical volume and G4PVPlacement for storing as the final
    //output
    G4ThreeVector groundPlaneTranslate(0,0,0.5*(dp_housingDimZ)
                                       + dp_eps + dp_groundPlaneDimZ*0.5);
    G4VPhysicalVolume * phys_groundPlane =
      new G4PVPlacement(0,groundPlaneTranslate,log_groundPlane,"GroundPlane", 
                        worldLogical,false,0,checkOverlaps);
      
    G4VisAttributes* groundPlaneVisAtt =
      new G4VisAttributes(G4Colour(0.0,1.0,1.0,0.5));
    groundPlaneVisAtt->SetVisibility(true);
    log_groundPlane->SetVisAttributes(groundPlaneVisAtt);  
    
    G4LatticePhysical* phys_groundPlaneLattice =
      new G4LatticePhysical(log_aluminumLattice,dp_polycryElScatMFP_Al,
                            dp_scDelta0_Al,dp_scTeff_Al,dp_scDn_Al,
                            dp_scTauQPTrap_Al);
    phys_groundPlaneLattice->SetMillerOrientation(1,0,0); 
    LM->RegisterLattice(phys_groundPlane,phys_groundPlaneLattice);

            
    //Set up the logical border surface
    G4CMPLogicalBorderSurface * border_siliconChip_groundPlane =
      new G4CMPLogicalBorderSurface("border_siliconChip_groundPlane",
                                    phys_siliconChip, phys_groundPlane,
                                    fSiAlInterface);
    G4CMPLogicalBorderSurface * border_groundPlane_siliconChip =
      new G4CMPLogicalBorderSurface("border_siliconChip_groundPlane",
                                    phys_groundPlane, phys_siliconChip,
                                    fSiAlInterface);
    G4CMPLogicalBorderSurface * border_world_groundPlane =
      new G4CMPLogicalBorderSurface("border_world_groundPlane", fWorldPhys,
                                    phys_groundPlane, fAlVacInterface);
    G4CMPLogicalBorderSurface * border_groundPlane_world =
      new G4CMPLogicalBorderSurface("border_groundPlane_world",
                                    phys_groundPlane, fWorldPhys,
                                    fAlVacInterface);

    //----------------------------------------------------------------------
    //Now set up the transmission line
    if (dp_useTransmissionLine) {

      //Since it's within the ground plane exactly:
      G4ThreeVector transmissionLineTranslate(0,0,0.0);
      QuasiparticleTransmissionLine * tLine =
        new QuasiparticleTransmissionLine(0,transmissionLineTranslate,
                                          "TransmissionLine",log_groundPlane,
                                          false,0,LM,fLogicalLatticeContainer,
                                          fBorderContainer,checkOverlaps);
      G4LogicalVolume * log_tLine = tLine->GetLogicalVolume();
      G4VPhysicalVolume * phys_tLine = tLine->GetPhysicalVolume();



      //Now, if we're using the chip and ground plane AND the transmission line
      //This gets a bit hairy, since the transmission line is composite of both
      //Nb and vacuum. So we'll access the list of physical objects present in
      //it and link those one-by-one to the silicon chip and to the world
      for( int iSubVol = 0; iSubVol < tLine->GetListOfAllFundamentalSubVolumes().size(); ++iSubVol){
        std::cout << "TLine sub volume names (to be used for boundaries): "
                  << std::get<1>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol])
                  << " with material "
                  << std::get<0>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol])
                  << std::endl;

        //Set the si/vac and si/Al interfaces to be symmetric in both dimensions
        std::string tempName1a = "border_siliconChip_"
          + std::get<1>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]);
        std::string tempName1b = "border_"
          + std::get<1>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol])
          + "_siliconChip";	  

        if (std::get<0>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Vacuum") != std::string::npos) {	  
          G4CMPLogicalBorderSurface * border1_siliconChip_transmissionLineEmpty=
            new G4CMPLogicalBorderSurface(tempName1a, phys_siliconChip,std::get<2>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fSiVacInterface);  
          G4CMPLogicalBorderSurface * border2_siliconChip_transmissionLineEmpty=
            new G4CMPLogicalBorderSurface(tempName1b,std::get<2>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), phys_siliconChip, fSiVacInterface);
        }
	
        if (std::get<0>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Aluminum") != std::string::npos) {
          G4CMPLogicalBorderSurface * border1_siliconChip_transmissionLineConductor = new G4CMPLogicalBorderSurface(tempName1a, phys_siliconChip, std::get<2>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fSiAlInterface);
          G4CMPLogicalBorderSurface * border2_siliconChip_transmissionLineConductor = new G4CMPLogicalBorderSurface(tempName1b, std::get<2>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), phys_siliconChip, fSiAlInterface);
        }


        //Set the world/vac and world/Al interfaces to be symmetric in both
        //directions
        std::string tempName2a = "border_world_"
          + std::get<1>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]);
        std::string tempName2b = "border_"
          + std::get<1>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol])
          + "_world";	  

        if (std::get<0>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Vacuum") != std::string::npos) {
          G4CMPLogicalBorderSurface * border1_world_transmissionLineEmpty = new G4CMPLogicalBorderSurface(tempName2a, fWorldPhys, std::get<2>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fVacVacInterface);
          G4CMPLogicalBorderSurface * border2_world_transmissionLineEmpty = new G4CMPLogicalBorderSurface(tempName2b, std::get<2>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fWorldPhys, fVacVacInterface);
        }
	
        if (std::get<0>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Aluminum") != std::string::npos) {
          G4CMPLogicalBorderSurface * border1_world_transmissionLineConductor = new G4CMPLogicalBorderSurface(tempName2a, fWorldPhys, std::get<2>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fAlVacInterface);
          G4CMPLogicalBorderSurface * border2_world_transmissionLineConductor = new G4CMPLogicalBorderSurface(tempName2b, std::get<2>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fWorldPhys, fAlVacInterface);
        }


        //Space here for linking to GROUND PLANE. Only need to link vacuum ones
        //since they're fully enclosing TL
        std::string tempName3a = "border_groundPlane_"
          + std::get<1>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]);
        std::string tempName3b = "border_"
          + std::get<1>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol])
          + "_groundPlane";
	
        if (std::get<0>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Vacuum") != std::string::npos) {
          G4CMPLogicalBorderSurface * border1_groundPlane_transmissionLineEmpty = new G4CMPLogicalBorderSurface(tempName3a, phys_groundPlane, std::get<2>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fAlVacInterface);
          G4CMPLogicalBorderSurface * border2_groundPlane_transmissionLineEmpty = new G4CMPLogicalBorderSurface(tempName3b, std::get<2>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), phys_groundPlane, fAlVacInterface);
        }	  
      }
    }


    //----------------------------------------------------------------------
    //Now set up a set of 6 resonator assemblies
    if (dp_useResonatorAssembly) {
      int nR = 6;
      for (int iR = 0; iR < nR; ++iR) {
      
        //First, get the translation vector for the resonator assembly
        //For the top three, don't do a rotation. For the bottom three, do
        G4ThreeVector resonatorAssemblyTranslate(0,0,0);
        G4RotationMatrix * rotAssembly = 0;
        if (iR <= 2) {
          resonatorAssemblyTranslate =
            G4ThreeVector(dp_resonatorLateralSpacing*(iR-1)
                          +dp_centralResonatorOffsetX,
                          0.5 * dp_resonatorAssemblyBaseAlDimY
                          + 0.5 * dp_transmissionLineCavityFullWidth,
                          0.0);
          rotAssembly = 0;
        }
        else {
          //Negative X offset because qubit is mirrored on underside
          resonatorAssemblyTranslate =
            G4ThreeVector(dp_resonatorLateralSpacing*(iR-4)
                          -dp_centralResonatorOffsetX, 
                          -1*(0.5 * dp_resonatorAssemblyBaseAlDimY
                              + 0.5 * dp_transmissionLineCavityFullWidth),
                          0.0);
          rotAssembly = new G4RotationMatrix();
          rotAssembly->rotateZ(180*deg);
        }
	
        char name[400];
        sprintf(name,"ResonatorAssembly_%d",iR);
        G4String resonatorAssemblyName(name);
        QuasiparticleResonatorAssembly * resonatorAssembly =
          new QuasiparticleResonatorAssembly(rotAssembly,
                                             resonatorAssemblyTranslate,
                                             resonatorAssemblyName,
                                             log_groundPlane,false,0,LM,
                                             fLogicalLatticeContainer,
                                             fBorderContainer,
                                             checkOverlaps);
        G4LogicalVolume * log_resonatorAssembly =
          resonatorAssembly->GetLogicalVolume();
        G4VPhysicalVolume * phys_resonatorAssembly =
          resonatorAssembly->GetPhysicalVolume();
	
        //Do the logical border creation now
        for( int iSubVol = 0; iSubVol < resonatorAssembly->GetListOfAllFundamentalSubVolumes().size(); ++iSubVol){
          std::cout << "TLine sub volume names (to be used for boundaries): " << std::get<1>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]) << " with material " << std::get<0>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]) << std::endl;
	  
	  
          //Set the chip/vacuum interfaces
          if (std::get<0>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Vacuum") != std::string::npos) {
            std::string tempName1 = "border_siliconChip_" + std::get<1>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]);
            std::string tempName2 = "border_" + std::get<1>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]) + "_siliconChip";
            G4CMPLogicalBorderSurface * border_siliconChip_resonatorAssemblyEmpty = new G4CMPLogicalBorderSurface(tempName1, std::get<2>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]), phys_siliconChip, fSiVacInterface);
            G4CMPLogicalBorderSurface * border_resonatorAssemblyEmpty_siliconChip = new G4CMPLogicalBorderSurface(tempName2, phys_siliconChip, std::get<2>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]), fSiVacInterface);
          }

          //Set the chip/aluminum interfaces
          if (std::get<0>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Aluminum") != std::string::npos) {
            std::string tempName1 = "border_siliconChip_" + std::get<1>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]);
            std::string tempName2 = "border_" + std::get<1>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]) + "_siliconChip";
            G4CMPLogicalBorderSurface * border_siliconChip_resonatorAssemblyConductor = new G4CMPLogicalBorderSurface(tempName1, phys_siliconChip, std::get<2>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]), fSiAlInterface);
            G4CMPLogicalBorderSurface * border_resonatorAssemblyConductor_siliconChip = new G4CMPLogicalBorderSurface(tempName2, std::get<2>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]), phys_siliconChip, fSiAlInterface);
          }

          //Set the world/vacuum interfaces (probably not necessary but
          //whatever, better to have everything well-defined
          if (std::get<0>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Vacuum") != std::string::npos) {
            std::string tempName1 = "border_world_" + std::get<1>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]);
            std::string tempName2 = "border_" + std::get<1>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]) + "_world";
            G4CMPLogicalBorderSurface * border_world_resonatorAssemblyEmpty = new G4CMPLogicalBorderSurface(tempName1, fWorldPhys, std::get<2>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]), fVacVacInterface);
            G4CMPLogicalBorderSurface * border_resonatorAssemblyEmpty_world = new G4CMPLogicalBorderSurface(tempName2, std::get<2>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]), fWorldPhys, fVacVacInterface);
          }

          //Set the world/aluminum interfaces
          if (std::get<0>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Aluminum") != std::string::npos) {
            std::string tempName1 = "border_world_" + std::get<1>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]);
            std::string tempName2 = "border_" + std::get<1>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]) + "_world";
            G4CMPLogicalBorderSurface * border_world_resonatorAssemblyConductor = new G4CMPLogicalBorderSurface(tempName1, fWorldPhys, std::get<2>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]), fAlVacInterface);
            G4CMPLogicalBorderSurface * border_resonatorAssemblyConductor_world = new G4CMPLogicalBorderSurface(tempName2, std::get<2>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]), fWorldPhys, fAlVacInterface);
          }

          //Set the TLcouplerConductor interface with the ground plane
          if (std::get<1>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("tlCouplingConductor") != std::string::npos) {
            std::string tempName1 = "border_groundPlane_" + std::get<1>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]);
            std::string tempName2 = "border_" + std::get<1>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]) + "_groundPlane";
            G4CMPLogicalBorderSurface * border_groundPlane_tlCouplingConductor = new G4CMPLogicalBorderSurface(tempName1, phys_groundPlane, std::get<2>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]), fAlAlInterface);
            G4CMPLogicalBorderSurface * border_tlCouplingConductor_groundPlane = new G4CMPLogicalBorderSurface(tempName2, std::get<2>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]), phys_groundPlane, fAlAlInterface);
          }


          //Set the TLCouplerEmpty interface with the ground plane
          if (std::get<1>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("tlCouplingEmpty") != std::string::npos) {
            std::string tempName1 = "border_groundPlane_" + std::get<1>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]);
            std::string tempName2 = "border_" + std::get<1>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]) + "_groundPlane";
            G4CMPLogicalBorderSurface * border_groundPlane_tlCouplingEmpty = new G4CMPLogicalBorderSurface(tempName1, phys_groundPlane, std::get<2>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]), fAlVacInterface);
            G4CMPLogicalBorderSurface * border_tlCouplingEmpty_groundPlane = new G4CMPLogicalBorderSurface(tempName2, std::get<2>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]), phys_groundPlane, fAlVacInterface);
          }

          //Set the baseLayer/groundplane interface
          if (std::get<1>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]) == resonatorAssemblyName) {
            std::string tempName1 = "border_groundPlane_" + std::get<1>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]);
            std::string tempName2 = "border_" + std::get<1>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]) + "_groundPlane";
            G4CMPLogicalBorderSurface * border_groundPlane_baselayer = new G4CMPLogicalBorderSurface(tempName1, phys_groundPlane, std::get<2>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]), fAlAlInterface);
            G4CMPLogicalBorderSurface * border_baselayer_groundPlane = new G4CMPLogicalBorderSurface(tempName2, std::get<2>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]), phys_groundPlane, fAlAlInterface);	      
          }
        }
      }
    }
  }
}





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Attach material properties and electrode/sensor handler to surface

void QuasiparticleDetectorConstruction::
AttachPhononSensor(G4CMPSurfaceProperty* surfProp) {
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


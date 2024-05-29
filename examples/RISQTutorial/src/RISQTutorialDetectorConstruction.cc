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

#include "RISQTutorialDetectorConstruction.hh"
#include "RISQTutorialSensitivity.hh"
#include "RISQTutorialQubitHousing.hh"
#include "RISQTutorialPad.hh"
#include "RISQTutorialTransmissionLine.hh"
#include "RISQTutorialStraightFluxLine.hh"
#include "RISQTutorialCornerFluxLine.hh"
#include "RISQTutorialResonatorAssembly.hh"
#include "G4CMPPhononElectrode.hh"
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

using namespace RISQTutorialDetectorParameters;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RISQTutorialDetectorConstruction::RISQTutorialDetectorConstruction()
  : fLiquidHelium(0), fGermanium(0), fAluminum(0), fTungsten(0),
    fWorldPhys(0), 
    fSuperconductorSensitivity(0), fConstructed(false) {;}//, fIfField(true) {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RISQTutorialDetectorConstruction::~RISQTutorialDetectorConstruction() {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* RISQTutorialDetectorConstruction::Construct()
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

void RISQTutorialDetectorConstruction::DefineMaterials()
{ 
  G4NistManager* nistManager = G4NistManager::Instance();

  fLiquidHelium = nistManager->FindOrBuildMaterial("G4_AIR"); // to be corrected
  fGermanium = nistManager->FindOrBuildMaterial("G4_Ge");
  fSilicon = nistManager->FindOrBuildMaterial("G4_Si");
  fAluminum = nistManager->FindOrBuildMaterial("G4_Al");
  fTungsten = nistManager->FindOrBuildMaterial("G4_W");
  fNiobium = nistManager->FindOrBuildMaterial("G4_Nb");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RISQTutorialDetectorConstruction::SetupGeometry()
{



  //---------------------------------------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------------------------------------
  // First, define border surface properties that can be referenced later
  const G4double GHz = 1e9 * hertz; 

  //the following coefficients and cutoff values are not well-motivated
  //the code below is used only to demonstrate how to set these values.
  const std::vector<G4double> anhCoeffs = {0,0,0,0,0,0};//Turn this off temporarily
  const std::vector<G4double> diffCoeffs = {1,0,0,0,0,0};//Explicitly make this 1 for now
  const std::vector<G4double> specCoeffs = {0,0,0,0,0,0};//Turn this off temporarily
  const G4double anhCutoff = 520., reflCutoff = 350.;   // Units external
    
  
  //These are just the definitions of the interface TYPES, not the interfaces themselves. These must be called in a set of loops
  //below, and invoke these surface definitions.
  if( !fConstructed ){
    fSiNbInterface = new G4CMPSurfaceProperty("SiNbInterface",
					      1.0, 0.0, 0.0, 0.0,
					      0.1, 1.0, 0.0, 0.0);
    fSiCopperInterface = new G4CMPSurfaceProperty("SiCopperInterface",
						  1.0, 0.0, 0.0, 0.0,
						  1.0, 0.0, 0.0, 0.0 );
    fSiVacuumInterface = new G4CMPSurfaceProperty("SiVacuumInterface",
						  0.0, 1.0, 0.0, 0.0,
						  0.0, 1.0, 0.0, 0.0 );
    

    fSiNbInterface->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					    diffCoeffs, specCoeffs, GHz, GHz, GHz);  
    fSiCopperInterface->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
						diffCoeffs, specCoeffs, GHz, GHz, GHz);  
    fSiVacuumInterface->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
						diffCoeffs, specCoeffs, GHz, GHz, GHz);

    //Add a phonon sensor to the interface properties here.
    AttachPhononSensor(fSiNbInterface);
  }
















  //---------------------------------------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------------------------------------
  // Now we start constructing the various components and their interfaces  
  //     
  // World
  //
  G4VSolid* solid_world = new G4Box("World",55.*cm,55.*cm,55.*cm);
  G4LogicalVolume* log_world = new G4LogicalVolume(solid_world,fLiquidHelium,"World");
  //  worldLogical->SetUserLimits(new G4UserLimits(10*mm, DBL_MAX, DBL_MAX, 0, 0));
  log_world->SetVisAttributes(G4VisAttributes::Invisible);
  fWorldPhys = new G4PVPlacement(0,
				 G4ThreeVector(),
				 log_world,
				 "World",
				 0,
                                 false,
				 0);
  
  
  bool checkOverlaps = true;

  



  //-------------------------------------------------------------------------------------------------------------------
  //First, set up the qubit chip substrate. By default, assume that we're using this. Otherwise, it's hard to establish
  //a sensitivity object for this.
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
							   log_world,
							   false,
							   0,
							   checkOverlaps);

  G4VisAttributes* siliconChipVisAtt= new G4VisAttributes(G4Colour(0.5,0.5,0.5));
  siliconChipVisAtt->SetVisibility(true);
  log_siliconChip->SetVisAttributes(siliconChipVisAtt);



  //Set up the G4CMP silicon lattice information using the G4LatticeManager
  // G4LatticeManager gives physics processes access to lattices by volume
  G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();
  G4LatticeLogical* log_siliconLattice = LM->LoadLattice(fSilicon, "Si");
    
  // G4LatticePhysical assigns G4LatticeLogical a physical orientation
  G4LatticePhysical* phys_siliconLattice = new G4LatticePhysical(log_siliconLattice);
  phys_siliconLattice->SetMillerOrientation(1,0,0); 
  LM->RegisterLattice(phys_siliconChip,phys_siliconLattice);

  //Set up border surfaces
  G4CMPLogicalBorderSurface * border_siliconChip_world = new G4CMPLogicalBorderSurface("border_siliconChip_world", phys_siliconChip, fWorldPhys, fSiVacuumInterface);

    



  //-------------------------------------------------------------------------------------------------------------------
  //If desired, set up the copper qubit housing
  if( dp_useQubitHousing ){
      
      
    RISQTutorialQubitHousing * qubitHousing = new RISQTutorialQubitHousing(0,
									   G4ThreeVector(0,0,0),
									   "QubitHousing",
									   log_world,
									   false,
									   0,
									   checkOverlaps);
    G4LogicalVolume * log_qubitHousing = qubitHousing->GetLogicalVolume();
    G4VPhysicalVolume * phys_qubitHousing = qubitHousing->GetPhysicalVolume();
    
    //Set up the logical border surface
    G4CMPLogicalBorderSurface * border_siliconChip_qubitHousing = new G4CMPLogicalBorderSurface("border_siliconChip_qubitHousing", phys_siliconChip, phys_qubitHousing, fSiCopperInterface);
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
							     log_world,
							     false,
							     0,
							     checkOverlaps);
    
    G4VisAttributes* groundPlaneVisAtt= new G4VisAttributes(G4Colour(0.0,1.0,1.0,0.5));
    groundPlaneVisAtt->SetVisibility(true);
    log_groundPlane->SetVisAttributes(groundPlaneVisAtt);
    
    
    //Set up the logical border surface
    G4CMPLogicalBorderSurface * border_siliconChip_groundPlane = new G4CMPLogicalBorderSurface("border_siliconChip_groundPlane", phys_siliconChip, phys_groundPlane, fSiNbInterface);


    

    //-------------------------------------------------------------------------------------------------------------------
    //Now set up the transmission line
    if( dp_useTransmissionLine ){
	
      G4ThreeVector transmissionLineTranslate(0,0,0.0);//Since it's within the ground plane exactly; 0.5*(dp_housingDimZ) + dp_eps + dp_groundPlaneDimZ*0.5 ); 
      RISQTutorialTransmissionLine * tLine = new RISQTutorialTransmissionLine(0,
									      transmissionLineTranslate,
									      "TransmissionLine",
									      log_groundPlane,
									      false,
									      0,
									      checkOverlaps);
      G4LogicalVolume * log_tLine = tLine->GetLogicalVolume();
      G4VPhysicalVolume * phys_tLine = tLine->GetPhysicalVolume();



      //Now, if we're using the chip and ground plane AND the transmission line
      //This gets a bit hairy, since the transmission line is composite of both Nb and vacuum.
      //So we'll access the list of physical objects present in it and link those one-by-one to the
      //silicon chip.
      for( int iSubVol = 0; iSubVol < tLine->GetListOfAllFundamentalSubVolumes().size(); ++iSubVol){
	std::cout << "TLine sub volume names (to be used for boundaries): " << std::get<1>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]) << " with material " << std::get<0>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]) << std::endl;

	std::string tempName = "border_siliconChip_" + std::get<1>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]);	  
	if( std::get<0>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Vacuum") != std::string::npos ){
	  G4CMPLogicalBorderSurface * border_siliconChip_transmissionLineEmpty = new G4CMPLogicalBorderSurface(tempName, phys_siliconChip, std::get<2>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fSiVacuumInterface);
	}
	if( std::get<0>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Niobium") != std::string::npos ){
	  G4CMPLogicalBorderSurface * border_siliconChip_transmissionLineConductor = new G4CMPLogicalBorderSurface(tempName, phys_siliconChip, std::get<2>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fSiNbInterface);
	}
      }
    }


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




  //---------------------------------------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------------------------------------
  // Now we establish a sensitivity object
  
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  if (!fSuperconductorSensitivity)
    fSuperconductorSensitivity = new RISQTutorialSensitivity("PhononElectrode");
  SDman->AddNewDetector(fSuperconductorSensitivity);
  log_siliconChip->SetSensitiveDetector(fSuperconductorSensitivity);



}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Set up a phonon sensor for this surface property object. I'm pretty sure that this
// phonon sensor doesn't get stapled to individual geometrical objects, but rather gets
// stapled to a surface property, but I'm not sure... have to ask mKelsey
void RISQTutorialDetectorConstruction::AttachPhononSensor(G4CMPSurfaceProperty * surfProp)
{
  //If no surface, don't do anything
  if(!surfProp) return;

  //Specify properties of the niobium sensors
  auto sensorProp = surfProp->GetPhononMaterialPropertiesTablePointer();
  sensorProp->AddConstProperty("filmAbsorption",0.0);              //NOT WELL MOTIVATED - probably parametrize and put on slider?
  sensorProp->AddConstProperty("filmThickness",90.*CLHEP::nm);     //Accurate for our thin film.
  sensorProp->AddConstProperty("gapEnergy",1.6e-3*CLHEP::eV);       //Reasonably motivated. Actually, looks like Novotny and Meincke are quoting 2Delta, and this is delta. Nuss and Goossen mention that Nb has a delta value closer to this.
  sensorProp->AddConstProperty("lowQPLimit",3.);                   //NOT WELL MOTIVATED YET -- Dunno how to inform this...
  sensorProp->AddConstProperty("phononLifetime",4.17*CLHEP::ps);   //Kaplan paper says 242ps for Al, same table says 4.17ps for characteristic time for Nb.
  sensorProp->AddConstProperty("phononLifetimeSlope",0.29);        //Based on guessing from Kaplan paper, I think this is material-agnostic?
  sensorProp->AddConstProperty("vSound",3.480*CLHEP::km/CLHEP::s); //True for room temperature, probably good to 10%ish - should follow up
  sensorProp->AddConstProperty("subgapAbsorption",0.0);            //Assuming that since we're mostly sensitive to quasiparticle density, phonon "heat" here isn't something that we're sensitive to? Unsure how to select this.

  //  sensorProp->AddConstProperty("gapEnergy",3.0e-3*CLHEP::eV);      //Reasonably motivated. Novotny and Meincke, 1975 (2.8-3.14 meV)
  //  sensorProp->AddConstProperty("phononLifetime",242.*ps);      //Kaplan paper says 242ps for Al, same table says 4.17ps for characteristic time for Nb.
  
  surfProp->SetPhononElectrode(new G4CMPPhononElectrode);
  
}

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
//

//Includes (basic)
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

//Includes (specific to this project)
#include "RISQTutorialResonatorAssembly.hh"
#include "RISQTutorialPad.hh"
#include "RISQTutorialDetectorParameters.hh"

using namespace RISQTutorialDetectorParameters;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Primary Constructor
RISQTutorialResonatorAssembly::RISQTutorialResonatorAssembly(G4RotationMatrix * pRot,
								       const G4ThreeVector & tLate,
								       const G4String & pName,
								       G4LogicalVolume * pMotherLogical,
								       G4bool pMany,
								       G4int pCopyNo,
								       G4bool pSurfChk)
{
  //Here, use the inputs to this to set up the geometry and fill out the PVPlacement data member,
  //which is the real output from this class (and which we'll access in our detector construction file
  //file.)

  ConstructResonatorAssembly(pRot,
			    tLate,
			    pName,
			    pMotherLogical,
			    pMany,
			    pCopyNo,
			    pSurfChk);

  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Default Constructor
RISQTutorialResonatorAssembly::RISQTutorialResonatorAssembly()
{  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Destructor
RISQTutorialResonatorAssembly::~RISQTutorialResonatorAssembly()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//Moving implementation down here so it's not in the constructor
void RISQTutorialResonatorAssembly::ConstructResonatorAssembly(G4RotationMatrix * pRot,
							     const G4ThreeVector & tLate,
							     const G4String & pName,
							     G4LogicalVolume * pMotherLogical,
							     G4bool pMany,
							     G4int pCopyNo,
							     G4bool pSurfChk)
{

  //Start with some preliminaries - NIST manager
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* niobium_mat = nist->FindOrBuildMaterial("G4_Nb");
  G4Material* air_mat = nist->FindOrBuildMaterial("G4_AIR");
  bool checkOverlaps = true;

  //Set up the niobium visualization
  G4VisAttributes* niobium_vis= new G4VisAttributes(G4Colour(0.0,1.0,1.0,0.5));
  niobium_vis->SetVisibility(true);
  G4VisAttributes* air_vis= new G4VisAttributes(G4Colour(0.5,0.5,0.5,0.5));
  air_vis->SetVisibility(true);
  




  

  //------------------------------------------------------------------------------------------
  //Start with a base layer of niobium into which our objects will fit. We'll return this in the end.
  G4String baseNbLayerName = pName;
  G4String baseNbLayerNameSolid = pName + "_solid";
  G4String baseNbLayerNameLog = pName + "_log";
  G4Box * solid_baseNbLayer = new G4Box(baseNbLayerNameSolid,
					0.5 * dp_resonatorAssemblyBaseNbDimX,
					0.5 * dp_resonatorAssemblyBaseNbDimY,
					0.5 * dp_resonatorAssemblyBaseNbDimZ);

  //Now attribute a physical material to the housing
  G4LogicalVolume * log_baseNbLayer = new G4LogicalVolume(solid_baseNbLayer,
							  niobium_mat,
							  baseNbLayerNameLog);
  log_baseNbLayer->SetVisAttributes(G4VisAttributes::Invisible);//niobium_vis);

  //Now, create a physical volume and G4PVPlacement for storing as the final output. This is the
  //top volume.
  G4VPhysicalVolume* phys_baseNbLayer = new G4PVPlacement(pRot,
							  tLate,
							  log_baseNbLayer,
							  baseNbLayerName, 
							  pMotherLogical,
							  pMany,
							  pCopyNo,
							  pSurfChk);

  //Also need an interface definition for this base layer...
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Niobium",baseNbLayerName,phys_baseNbLayer));



  //------------------------------------------------------------------------------------------
  //Now make the various components of the resonator array: line+coupling, shunt capacitance (cross), and qubit
  MakeResonatorLine(pName,log_baseNbLayer);
  MakeShuntCapacitorCross(pName,log_baseNbLayer);  

  /*  

  


  //------------------------------------------------------------------------------------------
  //Now make the pads of the transmission line. The logical mother volume of these is the base
  //Nb layer. Pads come with their own visualization attributes already set.
  
  //Pad 1:
  G4String pad1Name = pName + "_TransmissionLinePad1";
  RISQTutorialPad * pad1 = new RISQTutorialPad(0,
							       G4ThreeVector(dp_transmissionLinePad1Offset,0,0),
							       pad1Name,
							       log_baseNiLayer,
							       false,
							       0,
							       checkOverlaps);
  G4LogicalVolume * log_pad1 = pad1->GetLogicalVolume();
  G4VPhysicalVolume * phys_pad1 = pad1->GetPhysicalVolume();
  

  //Pad 2: rotate around Z axis by 180 degrees
  G4String pad2Name = pName + "_TransmissionLinePad2";
  G4RotationMatrix * pad2Rot = new G4RotationMatrix();
  pad2Rot->rotateZ(180*deg);
  RISQTutorialPad * pad2 = new RISQTutorialPad(pad2Rot,
							       G4ThreeVector(dp_transmissionLinePad2Offset,0,0),
							       pad2Name,
							       log_baseNiLayer,
							       false,
							       0,
							       checkOverlaps);
  G4LogicalVolume * log_pad2 = pad1->GetLogicalVolume();
  G4VPhysicalVolume * phys_pad2 = pad1->GetPhysicalVolume();

  


  
  //------------------------------------------------------------------------------------------
  //Now make the transmission line itself, in two parts: empty and conductor. The logical
  //mother volume of this is the base Ni layer.

  G4String tlNameEmpty = pName + "_TransmissionLineEmpty";
  G4String tlNameEmptySolid = tlNameEmpty + "_solid";
  G4String tlNameEmptyLog = tlNameEmpty + "_log";
  G4Box * solid_transmissionLineEmpty = new G4Box(tlNameEmptySolid,
						  0.5 * (dp_transmissionLinePad2Offset - dp_transmissionLinePad1Offset - 2 * dp_padEmptyPart2TrdZ - 2 * 0.5 * dp_padEmptyPart1DimX),
						  0.5 * dp_transmissionLineCavityFullWidth,
						  0.5 * dp_transmissionLineBaseLayerDimZ);

  G4LogicalVolume * log_transmissionLineEmpty = new G4LogicalVolume(solid_transmissionLineEmpty,
								    air_mat,
								    tlNameEmptyLog);
  log_transmissionLineEmpty->SetVisAttributes(air_vis);//G4VisAttributes::Invisible);
  
  G4VPhysicalVolume * phys_transmissionLineEmpty = new G4PVPlacement(0,
								     G4ThreeVector(0,0,0),
								     log_transmissionLineEmpty,
								     tlNameEmpty,
								     log_baseNiLayer,
								     false,
								     0,
								     true);


  
  G4String tlNameConductor = pName + "_TransmissionLineConductor";
  G4String tlNameConductorSolid = tlNameConductor + "_solid";
  G4String tlNameConductorLog = tlNameConductor + "_log";
  G4Box * solid_transmissionLineConductor = new G4Box(tlNameConductorSolid,
						      0.5 * (dp_transmissionLinePad2Offset - dp_transmissionLinePad1Offset - 2 * dp_padEmptyPart2TrdZ - 2 * 0.5 * dp_padEmptyPart1DimX),
						      0.5 * dp_transmissionLineConductorWidth,
						      0.5 * dp_transmissionLineBaseLayerDimZ);
  
  G4LogicalVolume * log_transmissionLineConductor = new G4LogicalVolume(solid_transmissionLineConductor,
									niobium_mat,
									tlNameConductorLog);
  log_transmissionLineConductor->SetVisAttributes(niobium_vis);
  
  G4VPhysicalVolume * phys_transmissionLineConductor = new G4PVPlacement(0,
									 G4ThreeVector(0,0,0),
									 log_transmissionLineConductor,
									 tlNameConductor,
									 log_transmissionLineEmpty,
									 false,
									 0,
									 true);

  


  */
  ///////////////////////////////////////////
  // Output logical/physical volume selection
  //-----------------------------------------

  fLog_output = log_baseNbLayer;
  fPhys_output = phys_baseNbLayer;


  /*

  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,0.647,0.0,0.9));
  simpleBoxVisAtt->SetVisibility(true);
  log_QubitHousing->SetVisAttributes(simpleBoxVisAtt);
  */


  

}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Make the resonator line
void RISQTutorialResonatorAssembly::MakeShuntCapacitorCross(G4String pName, G4LogicalVolume * log_baseNbLayer)
{


  //Materials and NIST
  //Start with some preliminaries - NIST manager
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* niobium_mat = nist->FindOrBuildMaterial("G4_Nb");
  G4Material* air_mat = nist->FindOrBuildMaterial("G4_AIR");
  bool checkOverlaps = true;
  
  //Set up the niobium visualization
  G4VisAttributes* niobium_vis= new G4VisAttributes(G4Colour(0.0,1.0,1.0,0.5));
  niobium_vis->SetVisibility(true);
  G4VisAttributes* air_vis= new G4VisAttributes(G4Colour(0.5,0.5,0.5,0.5));
  air_vis->SetVisibility(true);


  //This will be made in two batches: one for "empty" space and one for "conductor" space (the line itself)
  //Each batch has the following elements strung together, in order:
  //1. A vertical block
  //2. A horizontal block
  
  //Some useful translations relative to the center of the plane in which all of this is embedded
  G4ThreeVector brCornerOfBaseNbLayer(0.5*dp_resonatorAssemblyBaseNbDimX,-0.5*dp_resonatorAssemblyBaseNbDimY,0); //Bottom right corner relative to center of the plane
  
  //------------------------------------------------------
  //Vertical block (empty/cavity)
  G4String shuntEmptyName = pName + "_shuntEmpty";
  G4String shuntEmptyNameSolid = shuntEmptyName + "_solid";
  G4String shuntEmptyNameLog = shuntEmptyName + "_log";
  G4Box * solid_shuntVertBlockEmpty = new G4Box("shuntVertBlockEmpty",0.5 * dp_shuntVertBlockEmptyDimX,0.5 * dp_shuntVertBlockEmptyDimY,0.5 * dp_shuntVertBlockEmptyDimZ);
  G4Box * solid_shuntHorizontalBlockEmpty = new G4Box("ShuntHorizontalBlockEmpty",0.5 * dp_shuntHorizontalBlockEmptyDimX,0.5 * dp_shuntHorizontalBlockEmptyDimY,0.5 * dp_shuntHorizontalBlockEmptyDimZ);
  G4UnionSolid * solid_shuntEmpty = new G4UnionSolid(shuntEmptyNameSolid,solid_shuntVertBlockEmpty,solid_shuntHorizontalBlockEmpty,0,G4ThreeVector(0,0,0));

  
  G4LogicalVolume * log_shuntEmpty = new G4LogicalVolume(solid_shuntEmpty,air_mat,shuntEmptyNameLog);
  G4ThreeVector shuntWrtBRCorner(-1*dp_shuntCenterToBottomRightCornerOfBaseLayerDimX,dp_shuntCenterToBottomRightCornerOfBaseLayerDimY,0);
  G4VPhysicalVolume * shuntEmpty = new G4PVPlacement(0,shuntWrtBRCorner+brCornerOfBaseNbLayer,log_shuntEmpty,shuntEmptyName,log_baseNbLayer,false,0,true);
  log_shuntEmpty->SetVisAttributes(air_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Vacuum",shuntEmptyName,shuntEmpty));


  //------------------------------------------------------
  //Vertical block (conductor)
  G4String shuntConductorName = pName + "_shuntConductor";
  G4String shuntConductorNameSolid = shuntConductorName + "_solid";
  G4String shuntConductorNameLog = shuntConductorName + "_log";
  G4Box * solid_shuntVertBlockConductor = new G4Box("shuntVertBlockConductor",0.5 * dp_shuntVertBlockConductorDimX,0.5 * dp_shuntVertBlockConductorDimY,0.5 * dp_shuntVertBlockConductorDimZ);
  G4Box * solid_shuntHorizontalBlockConductor = new G4Box("ShuntHorizontalBlockConductor",0.5 * dp_shuntHorizontalBlockConductorDimX,0.5 * dp_shuntHorizontalBlockConductorDimY,0.5 * dp_shuntHorizontalBlockConductorDimZ);
  G4UnionSolid * solid_shuntConductor = new G4UnionSolid(shuntConductorNameSolid,solid_shuntVertBlockConductor,solid_shuntHorizontalBlockConductor,0,G4ThreeVector(0,0,0));  
  G4LogicalVolume * log_shuntConductor = new G4LogicalVolume(solid_shuntConductor,niobium_mat,shuntConductorNameLog);
  G4VPhysicalVolume * shuntConductor = new G4PVPlacement(0,G4ThreeVector(0,0,0),log_shuntConductor,shuntConductorName,log_shuntEmpty,false,0,true);
  log_shuntConductor->SetVisAttributes(niobium_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Niobium",shuntConductorName,shuntConductor));




}
  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Make the resonator line
void RISQTutorialResonatorAssembly::MakeResonatorLine(G4String pName, G4LogicalVolume * log_baseNbLayer)
{

  //Materials and NIST
  //Start with some preliminaries - NIST manager
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* niobium_mat = nist->FindOrBuildMaterial("G4_Nb");
  G4Material* air_mat = nist->FindOrBuildMaterial("G4_AIR");
  bool checkOverlaps = true;

  //Set up the niobium visualization
  G4VisAttributes* niobium_vis= new G4VisAttributes(G4Colour(0.0,1.0,1.0,0.5));
  niobium_vis->SetVisibility(true);
  G4VisAttributes* air_vis= new G4VisAttributes(G4Colour(0.5,0.5,0.5,0.5));
  air_vis->SetVisibility(true);

  
  
  //This will be made in two batches: one for "empty" space and one for "conductor" space (the line itself)
  //Each batch has the following elements strung together, in order:
  //1. Coupling to the transmission line - DONE
  //2. Curve 1 - DONE
  //3. Curve 2 - DONE
  //4. Straight horizontal line 1 - DONE
  //5. Half circle 1 - DONE
  //6. Straight horizontal line 2 - DONE
  //7. Half circle 2 - DONE
  //8. Straight horizontal line 3 - DONE
  //9. Half circle 3 - DONE
  //10. - DONE
  //11. Half circle 4 - DONE
  //12. - DONE
  //13. Half circle 5 - DONE
  //14. - DONE
  //15. Half circle 6 - DONE
  //15. Straight horizontal line 
  //16. Curve 3
  //17. Vertical straight to shunt coupler
  //18. Shunt coupler horizontal
  //19. Shunt coupler upper left lobe
  //20. Shunt coupler upper right lobe


  //Some useful translations relative to the center of the plane in which all of this is embedded
  G4ThreeVector brCornerOfBaseNbLayer(0.5*dp_resonatorAssemblyBaseNbDimX,-0.5*dp_resonatorAssemblyBaseNbDimY,0); //Bottom right corner relative to center of the plane
  
  

  //------------------------------------------------------
  //Coupling to the transmission line (empty/cavity)
  G4String tlCouplingEmptyName = pName + "_tlCouplingEmpty";
  G4String tlCouplingEmptyNameSolid = tlCouplingEmptyName + "_solid";
  G4String tlCouplingEmptyNameLog = tlCouplingEmptyName + "_log";
  G4Box * solid_tlCouplingEmpty = new G4Box(tlCouplingEmptyNameSolid,0.5 * dp_tlCouplingEmptyDimX,0.5 * dp_tlCouplingEmptyDimY,0.5 * dp_tlCouplingEmptyDimZ);
  G4LogicalVolume * log_tlCouplingEmpty = new G4LogicalVolume(solid_tlCouplingEmpty,air_mat,tlCouplingEmptyNameLog);
  G4ThreeVector tlCouplingWrtBRCorner(-0.5*dp_tlCouplingEmptyDimX,0.5*dp_tlCouplingEmptyDimY + dp_resonatorAssemblyBaseNbEdgeBottomDimY,0.0); //Good for empty or conductor
  G4VPhysicalVolume * tlCouplingEmpty = new G4PVPlacement(0,tlCouplingWrtBRCorner+brCornerOfBaseNbLayer,log_tlCouplingEmpty,tlCouplingEmptyName,log_baseNbLayer,false,0,true);
  log_tlCouplingEmpty->SetVisAttributes(air_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Vacuum",tlCouplingEmptyName,tlCouplingEmpty));


  
  //------------------------------------------------------
  //Coupling to the transmission line (conductor)
  G4String tlCouplingConductorName = pName + "_tlCouplingConductor";
  G4String tlCouplingConductorNameSolid = tlCouplingConductorName + "_solid";
  G4String tlCouplingConductorNameLog = tlCouplingConductorName + "_log";
  G4Box * solid_tlCouplingConductor = new G4Box(tlCouplingConductorNameSolid,0.5 * dp_tlCouplingConductorDimX,0.5 * dp_tlCouplingConductorDimY,0.5 * dp_tlCouplingConductorDimZ);
  G4LogicalVolume * log_tlCouplingConductor = new G4LogicalVolume(solid_tlCouplingConductor,niobium_mat,tlCouplingConductorNameLog);
  G4VPhysicalVolume * tlCouplingConductor = new G4PVPlacement(0,G4ThreeVector(0,0,0),log_tlCouplingConductor,tlCouplingConductorName,log_tlCouplingEmpty,false,0,true);
  log_tlCouplingConductor->SetVisAttributes(niobium_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Niobium",tlCouplingConductorName,tlCouplingConductor));


  //------------------------------------------------------
  //Curve 1 (empty/cavity)
  G4String curve1EmptyName = pName + "_curve1Empty";
  G4String curve1EmptyNameSolid = curve1EmptyName + "_solid";
  G4String curve1EmptyNameLog = curve1EmptyName + "_log";
  G4Tubs * solid_curve1Empty = new G4Tubs(curve1EmptyNameSolid,dp_resonatorAssemblyCurveSmallestRadius,dp_resonatorAssemblyCurveSmallestRadius + dp_tlCouplingEmptyDimY,dp_curveEmptyDimZ/2.0,180.*deg,90.*deg);
  G4LogicalVolume * log_curve1Empty = new G4LogicalVolume(solid_curve1Empty,air_mat,curve1EmptyNameLog);  
  G4ThreeVector curve1WrtBRCorner(-1*dp_tlCouplingEmptyDimX,0.5*dp_tlCouplingEmptyDimY + dp_resonatorAssemblyBaseNbEdgeBottomDimY + dp_resonatorAssemblyCurveCentralRadius,0.0); //Good for empty or conductor
  G4VPhysicalVolume * curve1Empty = new G4PVPlacement(0,curve1WrtBRCorner+brCornerOfBaseNbLayer,log_curve1Empty,curve1EmptyName,log_baseNbLayer,false,0,true);
  log_curve1Empty->SetVisAttributes(air_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Vacuum",curve1EmptyName,curve1Empty));

  //------------------------------------------------------
  //Curve 1 (conductor)
  G4String curve1ConductorName = pName + "_curve1Conductor";
  G4String curve1ConductorNameSolid = curve1ConductorName + "_solid";
  G4String curve1ConductorNameLog = curve1ConductorName + "_log";
  G4Tubs * solid_curve1Conductor = new G4Tubs(curve1ConductorNameSolid,dp_resonatorAssemblyCurveCentralRadius - dp_tlCouplingConductorDimY/2.0,dp_resonatorAssemblyCurveCentralRadius + dp_tlCouplingConductorDimY/2.0,dp_curveEmptyDimZ/2.0,180.*deg,90.*deg);
  G4LogicalVolume * log_curve1Conductor = new G4LogicalVolume(solid_curve1Conductor,niobium_mat,curve1ConductorNameLog);  
  G4VPhysicalVolume * curve1Conductor = new G4PVPlacement(0,G4ThreeVector(0,0,0),log_curve1Conductor,curve1ConductorName,log_curve1Empty,false,0,true);
  log_curve1Conductor->SetVisAttributes(niobium_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Niobium",curve1ConductorName,curve1Conductor));


  //------------------------------------------------------
  //Curve 2 (empty/cavity)
  G4String curve2EmptyName = pName + "_curve2Empty";
  G4String curve2EmptyNameSolid = curve2EmptyName + "_solid";
  G4String curve2EmptyNameLog = curve2EmptyName + "_log";
  G4Tubs * solid_curve2Empty = new G4Tubs(curve2EmptyNameSolid,dp_resonatorAssemblyCurveSmallestRadius,dp_resonatorAssemblyCurveSmallestRadius + dp_tlCouplingEmptyDimY,dp_curveEmptyDimZ/2.0,0.*deg,90.*deg);
  G4LogicalVolume * log_curve2Empty = new G4LogicalVolume(solid_curve2Empty,air_mat,curve2EmptyNameLog);  
  G4ThreeVector curve2WrtBRCorner(-1*dp_tlCouplingEmptyDimX - 2*dp_resonatorAssemblyCurveCentralRadius,0.5*dp_tlCouplingEmptyDimY + dp_resonatorAssemblyBaseNbEdgeBottomDimY + dp_resonatorAssemblyCurveCentralRadius,0.0); //Good for empty or conductor
  G4VPhysicalVolume * curve2Empty = new G4PVPlacement(0,curve2WrtBRCorner+brCornerOfBaseNbLayer,log_curve2Empty,curve2EmptyName,log_baseNbLayer,false,0,true);
  log_curve2Empty->SetVisAttributes(air_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Vacuum",curve2EmptyName,curve2Empty));

  

  //------------------------------------------------------
  //Curve 2 (conductor)
  G4String curve2ConductorName = pName + "_curve2Conductor";
  G4String curve2ConductorNameSolid = curve2ConductorName + "_solid";
  G4String curve2ConductorNameLog = curve2ConductorName + "_log";
  G4Tubs * solid_curve2Conductor = new G4Tubs(curve2ConductorNameSolid,dp_resonatorAssemblyCurveCentralRadius - dp_tlCouplingConductorDimY/2.0,dp_resonatorAssemblyCurveCentralRadius + dp_tlCouplingConductorDimY/2.0,dp_curveEmptyDimZ/2.0,0.*deg,90.*deg);
  G4LogicalVolume * log_curve2Conductor = new G4LogicalVolume(solid_curve2Conductor,niobium_mat,curve2ConductorNameLog);  
  G4VPhysicalVolume * curve2Conductor = new G4PVPlacement(0,G4ThreeVector(0,0,0),log_curve2Conductor,curve2ConductorName,log_curve2Empty,false,0,true);
  log_curve2Conductor->SetVisAttributes(niobium_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Niobium",curve2ConductorName,curve2Conductor));




  //------------------------------------------------------
  //Straight horizontal line (SHL) 1, (empty/cavity)
  G4String shl1EmptyName = pName + "_shl1Empty";
  G4String shl1EmptyNameSolid = shl1EmptyName + "_solid";
  G4String shl1EmptyNameLog = shl1EmptyName + "_log";
  G4Box * solid_shl1Empty = new G4Box(shl1EmptyNameSolid,0.5 * dp_shl1EmptyDimX,0.5 * dp_shl1EmptyDimY,0.5 * dp_shl1EmptyDimZ);
  G4LogicalVolume * log_shl1Empty = new G4LogicalVolume(solid_shl1Empty,air_mat,shl1EmptyNameLog);
  G4ThreeVector shl1WrtBRCorner = curve2WrtBRCorner + G4ThreeVector(-0.5*dp_shl1EmptyDimX,dp_resonatorAssemblyCurveCentralRadius,0);  
  G4VPhysicalVolume * shl1Empty = new G4PVPlacement(0,shl1WrtBRCorner+brCornerOfBaseNbLayer,log_shl1Empty,shl1EmptyName,log_baseNbLayer,false,0,true);
  log_shl1Empty->SetVisAttributes(air_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Vacuum",shl1EmptyName,shl1Empty));


  //------------------------------------------------------
  //Straight horizontal line (SHL1) 1 (conductor)
  G4String shl1ConductorName = pName + "_shl1Conductor";
  G4String shl1ConductorNameSolid = shl1ConductorName + "_solid";
  G4String shl1ConductorNameLog = shl1ConductorName + "_log";
  G4Box * solid_shl1Conductor = new G4Box(shl1ConductorNameSolid,0.5 * dp_shl1ConductorDimX,0.5 * dp_shl1ConductorDimY,0.5 * dp_shl1ConductorDimZ);
  G4LogicalVolume * log_shl1Conductor = new G4LogicalVolume(solid_shl1Conductor,niobium_mat,shl1ConductorNameLog);
  G4VPhysicalVolume * shl1Conductor = new G4PVPlacement(0,G4ThreeVector(0,0,0),log_shl1Conductor,shl1ConductorName,log_shl1Empty,false,0,true);
  log_shl1Conductor->SetVisAttributes(niobium_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Niobium",shl1ConductorName,shl1Conductor));




  //------------------------------------------------------
  //HalfCircle 1 (empty/cavity)
  G4String halfCircle1EmptyName = pName + "_halfCircle1Empty";
  G4String halfCircle1EmptyNameSolid = halfCircle1EmptyName + "_solid";
  G4String halfCircle1EmptyNameLog = halfCircle1EmptyName + "_log";
  G4Tubs * solid_halfCircle1Empty = new G4Tubs(halfCircle1EmptyNameSolid,dp_resonatorAssemblyCurveSmallestRadius,dp_resonatorAssemblyCurveSmallestRadius + dp_tlCouplingEmptyDimY,dp_curveEmptyDimZ/2.0,90.*deg,180.*deg);
  G4LogicalVolume * log_halfCircle1Empty = new G4LogicalVolume(solid_halfCircle1Empty,air_mat,halfCircle1EmptyNameLog);  
  G4ThreeVector halfCircle1WrtBRCorner = shl1WrtBRCorner + G4ThreeVector(-0.5*dp_shl1EmptyDimX,dp_resonatorAssemblyCurveCentralRadius,0.0);
  G4VPhysicalVolume * halfCircle1Empty = new G4PVPlacement(0,halfCircle1WrtBRCorner+brCornerOfBaseNbLayer,log_halfCircle1Empty,halfCircle1EmptyName,log_baseNbLayer,false,0,true);
  log_halfCircle1Empty->SetVisAttributes(air_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Vacuum",halfCircle1EmptyName,halfCircle1Empty));
  
  //------------------------------------------------------
  //HalfCircle 1 (conductor)
  G4String halfCircle1ConductorName = pName + "_halfCircle1Conductor";
  G4String halfCircle1ConductorNameSolid = halfCircle1ConductorName + "_solid";
  G4String halfCircle1ConductorNameLog = halfCircle1ConductorName + "_log";
  G4Tubs * solid_halfCircle1Conductor = new G4Tubs(halfCircle1ConductorNameSolid,dp_resonatorAssemblyCurveCentralRadius - dp_tlCouplingConductorDimY/2.0,dp_resonatorAssemblyCurveCentralRadius + dp_tlCouplingConductorDimY/2.0,dp_curveEmptyDimZ/2.0,90.*deg,180.*deg);
  G4LogicalVolume * log_halfCircle1Conductor = new G4LogicalVolume(solid_halfCircle1Conductor,niobium_mat,halfCircle1ConductorNameLog);  
  G4VPhysicalVolume * halfCircle1Conductor = new G4PVPlacement(0,G4ThreeVector(0,0,0),log_halfCircle1Conductor,halfCircle1ConductorName,log_halfCircle1Empty,false,0,true);
  log_halfCircle1Conductor->SetVisAttributes(niobium_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Niobium",halfCircle1ConductorName,halfCircle1Conductor));


  
  //------------------------------------------------------
  //Straight horizontal line (SHL) 2, (empty/cavity)
  G4String shl2EmptyName = pName + "_shl2Empty";
  G4String shl2EmptyNameSolid = shl2EmptyName + "_solid";
  G4String shl2EmptyNameLog = shl2EmptyName + "_log";
  G4Box * solid_shl2Empty = new G4Box(shl2EmptyNameSolid,0.5 * dp_shl2EmptyDimX,0.5 * dp_shl2EmptyDimY,0.5 * dp_shl2EmptyDimZ);
  G4LogicalVolume * log_shl2Empty = new G4LogicalVolume(solid_shl2Empty,air_mat,shl2EmptyNameLog);
  G4ThreeVector shl2WrtBRCorner = halfCircle1WrtBRCorner + G4ThreeVector(0.5*dp_shl2EmptyDimX,dp_resonatorAssemblyCurveCentralRadius,0);
  G4VPhysicalVolume * shl2Empty = new G4PVPlacement(0,shl2WrtBRCorner+brCornerOfBaseNbLayer,log_shl2Empty,shl2EmptyName,log_baseNbLayer,false,0,true);
  log_shl2Empty->SetVisAttributes(air_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Vacuum",shl2EmptyName,shl2Empty));

  //------------------------------------------------------
  //Straight horizontal line (SHL) 2 (conductor)
  G4String shl2ConductorName = pName + "_shl2Conductor";
  G4String shl2ConductorNameSolid = shl2ConductorName + "_solid";
  G4String shl2ConductorNameLog = shl2ConductorName + "_log";
  G4Box * solid_shl2Conductor = new G4Box(shl2ConductorNameSolid,0.5 * dp_shl2ConductorDimX,0.5 * dp_shl2ConductorDimY,0.5 * dp_shl2ConductorDimZ);
  G4LogicalVolume * log_shl2Conductor = new G4LogicalVolume(solid_shl2Conductor,niobium_mat,shl2ConductorNameLog);
  G4VPhysicalVolume * shl2Conductor = new G4PVPlacement(0,G4ThreeVector(0,0,0),log_shl2Conductor,shl2ConductorName,log_shl2Empty,false,0,true);
  log_shl2Conductor->SetVisAttributes(niobium_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Niobium",shl2ConductorName,shl2Conductor));


  
  //------------------------------------------------------
  //HalfCircle 2 (empty/cavity)
  G4String halfCircle2EmptyName = pName + "_halfCircle2Empty";
  G4String halfCircle2EmptyNameSolid = halfCircle2EmptyName + "_solid";
  G4String halfCircle2EmptyNameLog = halfCircle2EmptyName + "_log";
  G4Tubs * solid_halfCircle2Empty = new G4Tubs(halfCircle2EmptyNameSolid,dp_resonatorAssemblyCurveSmallestRadius,dp_resonatorAssemblyCurveSmallestRadius + dp_tlCouplingEmptyDimY,dp_curveEmptyDimZ/2.0,270.*deg,180.*deg);
  G4LogicalVolume * log_halfCircle2Empty = new G4LogicalVolume(solid_halfCircle2Empty,air_mat,halfCircle2EmptyNameLog);  
  G4ThreeVector halfCircle2WrtBRCorner = shl2WrtBRCorner + G4ThreeVector(0.5*dp_shl2EmptyDimX,dp_resonatorAssemblyCurveCentralRadius,0.0);
  G4VPhysicalVolume * halfCircle2Empty = new G4PVPlacement(0,halfCircle2WrtBRCorner+brCornerOfBaseNbLayer,log_halfCircle2Empty,halfCircle2EmptyName,log_baseNbLayer,false,0,true);
  log_halfCircle2Empty->SetVisAttributes(air_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Vacuum",halfCircle2EmptyName,halfCircle2Empty));

  //------------------------------------------------------
  //HalfCircle 2 (conductor)
  G4String halfCircle2ConductorName = pName + "_halfCircle2Conductor";
  G4String halfCircle2ConductorNameSolid = halfCircle2ConductorName + "_solid";
  G4String halfCircle2ConductorNameLog = halfCircle2ConductorName + "_log";
  G4Tubs * solid_halfCircle2Conductor = new G4Tubs(halfCircle2ConductorNameSolid,dp_resonatorAssemblyCurveCentralRadius - dp_tlCouplingConductorDimY/2.0,dp_resonatorAssemblyCurveCentralRadius + dp_tlCouplingConductorDimY/2.0,dp_curveEmptyDimZ/2.0,270.*deg,180.*deg);
  G4LogicalVolume * log_halfCircle2Conductor = new G4LogicalVolume(solid_halfCircle2Conductor,niobium_mat,halfCircle2ConductorNameLog);  
  G4VPhysicalVolume * halfCircle2Conductor = new G4PVPlacement(0,G4ThreeVector(0,0,0),log_halfCircle2Conductor,halfCircle2ConductorName,log_halfCircle2Empty,false,0,true);
  log_halfCircle2Conductor->SetVisAttributes(niobium_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Niobium",halfCircle2ConductorName,halfCircle2Conductor));
  
  //------------------------------------------------------
  //Straight horizontal line (SHL) 3, (empty/cavity)
  G4String shl3EmptyName = pName + "_shl3Empty";
  G4String shl3EmptyNameSolid = shl3EmptyName + "_solid";
  G4String shl3EmptyNameLog = shl3EmptyName + "_log";
  G4Box * solid_shl3Empty = new G4Box(shl3EmptyNameSolid,0.5 * dp_shl3EmptyDimX,0.5 * dp_shl3EmptyDimY,0.5 * dp_shl3EmptyDimZ);
  G4LogicalVolume * log_shl3Empty = new G4LogicalVolume(solid_shl3Empty,air_mat,shl3EmptyNameLog);
  G4ThreeVector shl3WrtBRCorner = halfCircle2WrtBRCorner + G4ThreeVector(-0.5*dp_shl3EmptyDimX,dp_resonatorAssemblyCurveCentralRadius,0);
  G4VPhysicalVolume * shl3Empty = new G4PVPlacement(0,shl3WrtBRCorner+brCornerOfBaseNbLayer,log_shl3Empty,shl3EmptyName,log_baseNbLayer,false,0,true);
  log_shl3Empty->SetVisAttributes(air_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Vacuum",shl3EmptyName,shl3Empty));

  //------------------------------------------------------
  //Straight horizontal line (SHL) 3 (conductor)
  G4String shl3ConductorName = pName + "_shl3Conductor";
  G4String shl3ConductorNameSolid = shl3ConductorName + "_solid";
  G4String shl3ConductorNameLog = shl3ConductorName + "_log";
  G4Box * solid_shl3Conductor = new G4Box(shl3ConductorNameSolid,0.5 * dp_shl3ConductorDimX,0.5 * dp_shl3ConductorDimY,0.5 * dp_shl3ConductorDimZ);
  G4LogicalVolume * log_shl3Conductor = new G4LogicalVolume(solid_shl3Conductor,niobium_mat,shl3ConductorNameLog);
  G4VPhysicalVolume * shl3Conductor = new G4PVPlacement(0,G4ThreeVector(0,0,0),log_shl3Conductor,shl3ConductorName,log_shl3Empty,false,0,true);
  log_shl3Conductor->SetVisAttributes(niobium_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Niobium",shl3ConductorName,shl3Conductor));

  //------------------------------------------------------
  //HalfCircle 3 (empty/cavity)
  G4String halfCircle3EmptyName = pName + "_halfCircle3Empty";
  G4String halfCircle3EmptyNameSolid = halfCircle3EmptyName + "_solid";
  G4String halfCircle3EmptyNameLog = halfCircle3EmptyName + "_log";
  G4Tubs * solid_halfCircle3Empty = new G4Tubs(halfCircle3EmptyNameSolid,dp_resonatorAssemblyCurveSmallestRadius,dp_resonatorAssemblyCurveSmallestRadius + dp_tlCouplingEmptyDimY,dp_curveEmptyDimZ/2.0,90.*deg,180.*deg);
  G4LogicalVolume * log_halfCircle3Empty = new G4LogicalVolume(solid_halfCircle3Empty,air_mat,halfCircle3EmptyNameLog);  
  G4ThreeVector halfCircle3WrtBRCorner = shl3WrtBRCorner + G4ThreeVector(-0.5*dp_shl3EmptyDimX,dp_resonatorAssemblyCurveCentralRadius,0.0);
  G4VPhysicalVolume * halfCircle3Empty = new G4PVPlacement(0,halfCircle3WrtBRCorner+brCornerOfBaseNbLayer,log_halfCircle3Empty,halfCircle3EmptyName,log_baseNbLayer,false,0,true);
  log_halfCircle3Empty->SetVisAttributes(air_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Vacuum",halfCircle3EmptyName,halfCircle3Empty));

  
  //------------------------------------------------------
  //HalfCircle 3 (conductor)
  G4String halfCircle3ConductorName = pName + "_halfCircle3Conductor";
  G4String halfCircle3ConductorNameSolid = halfCircle3ConductorName + "_solid";
  G4String halfCircle3ConductorNameLog = halfCircle3ConductorName + "_log";
  G4Tubs * solid_halfCircle3Conductor = new G4Tubs(halfCircle3ConductorNameSolid,dp_resonatorAssemblyCurveCentralRadius - dp_tlCouplingConductorDimY/2.0,dp_resonatorAssemblyCurveCentralRadius + dp_tlCouplingConductorDimY/2.0,dp_curveEmptyDimZ/2.0,90.*deg,180.*deg);
  G4LogicalVolume * log_halfCircle3Conductor = new G4LogicalVolume(solid_halfCircle3Conductor,niobium_mat,halfCircle3ConductorNameLog);  
  G4VPhysicalVolume * halfCircle3Conductor = new G4PVPlacement(0,G4ThreeVector(0,0,0),log_halfCircle3Conductor,halfCircle3ConductorName,log_halfCircle3Empty,false,0,true);
  log_halfCircle3Conductor->SetVisAttributes(niobium_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Niobium",halfCircle3ConductorName,halfCircle3Conductor));


  //------------------------------------------------------
  //Straight horizontal line (SHL) 4, (empty/cavity)
  G4String shl4EmptyName = pName + "_shl4Empty";
  G4String shl4EmptyNameSolid = shl4EmptyName + "_solid";
  G4String shl4EmptyNameLog = shl4EmptyName + "_log";
  G4Box * solid_shl4Empty = new G4Box(shl4EmptyNameSolid,0.5 * dp_shl4EmptyDimX,0.5 * dp_shl4EmptyDimY,0.5 * dp_shl4EmptyDimZ);
  G4LogicalVolume * log_shl4Empty = new G4LogicalVolume(solid_shl4Empty,air_mat,shl4EmptyNameLog);
  G4ThreeVector shl4WrtBRCorner = halfCircle3WrtBRCorner + G4ThreeVector(0.5*dp_shl4EmptyDimX,dp_resonatorAssemblyCurveCentralRadius,0);
  G4VPhysicalVolume * shl4Empty = new G4PVPlacement(0,shl4WrtBRCorner+brCornerOfBaseNbLayer,log_shl4Empty,shl4EmptyName,log_baseNbLayer,false,0,true);
  log_shl4Empty->SetVisAttributes(air_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Vacuum",shl4EmptyName,shl4Empty));

  //------------------------------------------------------
  //Straight horizontal line (SHL) 4 (conductor)
  G4String shl4ConductorName = pName + "_shl4Conductor";
  G4String shl4ConductorNameSolid = shl4ConductorName + "_solid";
  G4String shl4ConductorNameLog = shl4ConductorName + "_log";
  G4Box * solid_shl4Conductor = new G4Box(shl4ConductorNameSolid,0.5 * dp_shl4ConductorDimX,0.5 * dp_shl4ConductorDimY,0.5 * dp_shl4ConductorDimZ);
  G4LogicalVolume * log_shl4Conductor = new G4LogicalVolume(solid_shl4Conductor,niobium_mat,shl4ConductorNameLog);
  G4VPhysicalVolume * shl4Conductor = new G4PVPlacement(0,G4ThreeVector(0,0,0),log_shl4Conductor,shl4ConductorName,log_shl4Empty,false,0,true);
  log_shl4Conductor->SetVisAttributes(niobium_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Niobium",shl4ConductorName,shl4Conductor));

    //------------------------------------------------------
  //HalfCircle 4 (empty/cavity)
  G4String halfCircle4EmptyName = pName + "_halfCircle4Empty";
  G4String halfCircle4EmptyNameSolid = halfCircle4EmptyName + "_solid";
  G4String halfCircle4EmptyNameLog = halfCircle4EmptyName + "_log";
  G4Tubs * solid_halfCircle4Empty = new G4Tubs(halfCircle4EmptyNameSolid,dp_resonatorAssemblyCurveSmallestRadius,dp_resonatorAssemblyCurveSmallestRadius + dp_tlCouplingEmptyDimY,dp_curveEmptyDimZ/2.0,270.*deg,180.*deg);
  G4LogicalVolume * log_halfCircle4Empty = new G4LogicalVolume(solid_halfCircle4Empty,air_mat,halfCircle4EmptyNameLog);  
  G4ThreeVector halfCircle4WrtBRCorner = shl4WrtBRCorner + G4ThreeVector(+0.5*dp_shl4EmptyDimX,dp_resonatorAssemblyCurveCentralRadius,0.0);
  G4VPhysicalVolume * halfCircle4Empty = new G4PVPlacement(0,halfCircle4WrtBRCorner+brCornerOfBaseNbLayer,log_halfCircle4Empty,halfCircle4EmptyName,log_baseNbLayer,false,0,true);
  log_halfCircle4Empty->SetVisAttributes(air_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Vacuum",halfCircle4EmptyName,halfCircle4Empty));

  //------------------------------------------------------
  //HalfCircle 4 (conductor)
  G4String halfCircle4ConductorName = pName + "_halfCircle4Conductor";
  G4String halfCircle4ConductorNameSolid = halfCircle4ConductorName + "_solid";
  G4String halfCircle4ConductorNameLog = halfCircle4ConductorName + "_log";
  G4Tubs * solid_halfCircle4Conductor = new G4Tubs(halfCircle4ConductorNameSolid,dp_resonatorAssemblyCurveCentralRadius - dp_tlCouplingConductorDimY/2.0,dp_resonatorAssemblyCurveCentralRadius + dp_tlCouplingConductorDimY/2.0,dp_curveEmptyDimZ/2.0,270.*deg,180.*deg);
  G4LogicalVolume * log_halfCircle4Conductor = new G4LogicalVolume(solid_halfCircle4Conductor,niobium_mat,halfCircle4ConductorNameLog);  
  G4VPhysicalVolume * halfCircle4Conductor = new G4PVPlacement(0,G4ThreeVector(0,0,0),log_halfCircle4Conductor,halfCircle4ConductorName,log_halfCircle4Empty,false,0,true);
  log_halfCircle4Conductor->SetVisAttributes(niobium_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Niobium",halfCircle4ConductorName,halfCircle4Conductor));





  //------------------------------------------------------
  //Straight horizontal line (SHL) 5, (empty/cavity)
  G4String shl5EmptyName = pName + "_shl5Empty";
  G4String shl5EmptyNameSolid = shl5EmptyName + "_solid";
  G4String shl5EmptyNameLog = shl5EmptyName + "_log";
  G4Box * solid_shl5Empty = new G4Box(shl5EmptyNameSolid,0.5 * dp_shl5EmptyDimX,0.5 * dp_shl5EmptyDimY,0.5 * dp_shl5EmptyDimZ);
  G4LogicalVolume * log_shl5Empty = new G4LogicalVolume(solid_shl5Empty,air_mat,shl5EmptyNameLog);
  G4ThreeVector shl5WrtBRCorner = halfCircle4WrtBRCorner + G4ThreeVector(-0.5*dp_shl5EmptyDimX,dp_resonatorAssemblyCurveCentralRadius,0);
  G4VPhysicalVolume * shl5Empty = new G4PVPlacement(0,shl5WrtBRCorner+brCornerOfBaseNbLayer,log_shl5Empty,shl5EmptyName,log_baseNbLayer,false,0,true);
  log_shl5Empty->SetVisAttributes(air_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Vacuum",shl5EmptyName,shl5Empty));
 

  //------------------------------------------------------
  //Straight horizontal line (SHL) 5 (conductor)
  G4String shl5ConductorName = pName + "_shl5Conductor";
  G4String shl5ConductorNameSolid = shl5ConductorName + "_solid";
  G4String shl5ConductorNameLog = shl5ConductorName + "_log";
  G4Box * solid_shl5Conductor = new G4Box(shl5ConductorNameSolid,0.5 * dp_shl5ConductorDimX,0.5 * dp_shl5ConductorDimY,0.5 * dp_shl5ConductorDimZ);
  G4LogicalVolume * log_shl5Conductor = new G4LogicalVolume(solid_shl5Conductor,niobium_mat,shl5ConductorNameLog);
  G4VPhysicalVolume * shl5Conductor = new G4PVPlacement(0,G4ThreeVector(0,0,0),log_shl5Conductor,shl5ConductorName,log_shl5Empty,false,0,true);
  log_shl5Conductor->SetVisAttributes(niobium_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Niobium",shl5ConductorName,shl5Conductor));

  //------------------------------------------------------
  //HalfCircle 5 (empty/cavity)
  G4String halfCircle5EmptyName = pName + "_halfCircle5Empty";
  G4String halfCircle5EmptyNameSolid = halfCircle5EmptyName + "_solid";
  G4String halfCircle5EmptyNameLog = halfCircle5EmptyName + "_log";
  G4Tubs * solid_halfCircle5Empty = new G4Tubs(halfCircle5EmptyNameSolid,dp_resonatorAssemblyCurveSmallestRadius,dp_resonatorAssemblyCurveSmallestRadius + dp_tlCouplingEmptyDimY,dp_curveEmptyDimZ/2.0,90.*deg,180.*deg);
  G4LogicalVolume * log_halfCircle5Empty = new G4LogicalVolume(solid_halfCircle5Empty,air_mat,halfCircle5EmptyNameLog);  
  G4ThreeVector halfCircle5WrtBRCorner = shl5WrtBRCorner + G4ThreeVector(-0.5*dp_shl5EmptyDimX,dp_resonatorAssemblyCurveCentralRadius,0.0);
  G4VPhysicalVolume * halfCircle5Empty = new G4PVPlacement(0,halfCircle5WrtBRCorner+brCornerOfBaseNbLayer,log_halfCircle5Empty,halfCircle5EmptyName,log_baseNbLayer,false,0,true);
  log_halfCircle5Empty->SetVisAttributes(air_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Vacuum",halfCircle5EmptyName,halfCircle5Empty));

  //------------------------------------------------------
  //HalfCircle 5 (conductor)
  G4String halfCircle5ConductorName = pName + "_halfCircle5Conductor";
  G4String halfCircle5ConductorNameSolid = halfCircle5ConductorName + "_solid";
  G4String halfCircle5ConductorNameLog = halfCircle5ConductorName + "_log";
  G4Tubs * solid_halfCircle5Conductor = new G4Tubs(halfCircle5ConductorNameSolid,dp_resonatorAssemblyCurveCentralRadius - dp_tlCouplingConductorDimY/2.0,dp_resonatorAssemblyCurveCentralRadius + dp_tlCouplingConductorDimY/2.0,dp_curveEmptyDimZ/2.0,90.*deg,180.*deg);
  G4LogicalVolume * log_halfCircle5Conductor = new G4LogicalVolume(solid_halfCircle5Conductor,niobium_mat,halfCircle5ConductorNameLog);  
  G4VPhysicalVolume * halfCircle5Conductor = new G4PVPlacement(0,G4ThreeVector(0,0,0),log_halfCircle5Conductor,halfCircle5ConductorName,log_halfCircle5Empty,false,0,true);
  log_halfCircle5Conductor->SetVisAttributes(niobium_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Niobium",halfCircle5ConductorName,halfCircle5Conductor));









  //------------------------------------------------------
  //Straight horizontal line (SHL) 6, (empty/cavity)
  G4String shl6EmptyName = pName + "_shl6Empty";
  G4String shl6EmptyNameSolid = shl6EmptyName + "_solid";
  G4String shl6EmptyNameLog = shl6EmptyName + "_log";
  G4Box * solid_shl6Empty = new G4Box(shl6EmptyNameSolid,0.5 * dp_shl6EmptyDimX,0.5 * dp_shl6EmptyDimY,0.5 * dp_shl6EmptyDimZ);
  G4LogicalVolume * log_shl6Empty = new G4LogicalVolume(solid_shl6Empty,air_mat,shl6EmptyNameLog);
  G4ThreeVector shl6WrtBRCorner = halfCircle5WrtBRCorner + G4ThreeVector(0.5*dp_shl6EmptyDimX,dp_resonatorAssemblyCurveCentralRadius,0);
  G4VPhysicalVolume * shl6Empty = new G4PVPlacement(0,shl6WrtBRCorner+brCornerOfBaseNbLayer,log_shl6Empty,shl6EmptyName,log_baseNbLayer,false,0,true);
  log_shl6Empty->SetVisAttributes(air_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Vacuum",shl6EmptyName,shl6Empty));


  //------------------------------------------------------
  //Straight horizontal line (SHL) 6 (conductor)
  G4String shl6ConductorName = pName + "_shl6Conductor";
  G4String shl6ConductorNameSolid = shl6ConductorName + "_solid";
  G4String shl6ConductorNameLog = shl6ConductorName + "_log";
  G4Box * solid_shl6Conductor = new G4Box(shl6ConductorNameSolid,0.5 * dp_shl6ConductorDimX,0.5 * dp_shl6ConductorDimY,0.5 * dp_shl6ConductorDimZ);
  G4LogicalVolume * log_shl6Conductor = new G4LogicalVolume(solid_shl6Conductor,niobium_mat,shl6ConductorNameLog);
  G4VPhysicalVolume * shl6Conductor = new G4PVPlacement(0,G4ThreeVector(0,0,0),log_shl6Conductor,shl6ConductorName,log_shl6Empty,false,0,true);
  log_shl6Conductor->SetVisAttributes(niobium_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Niobium",shl6ConductorName,shl6Conductor));


  //------------------------------------------------------
  //HalfCircle 6 (empty/cavity)
  G4String halfCircle6EmptyName = pName + "_halfCircle6Empty";
  G4String halfCircle6EmptyNameSolid = halfCircle6EmptyName + "_solid";
  G4String halfCircle6EmptyNameLog = halfCircle6EmptyName + "_log";
  G4Tubs * solid_halfCircle6Empty = new G4Tubs(halfCircle6EmptyNameSolid,dp_resonatorAssemblyCurveSmallestRadius,dp_resonatorAssemblyCurveSmallestRadius + dp_tlCouplingEmptyDimY,dp_curveEmptyDimZ/2.0,270.*deg,180.*deg);
  G4LogicalVolume * log_halfCircle6Empty = new G4LogicalVolume(solid_halfCircle6Empty,air_mat,halfCircle6EmptyNameLog);  
  G4ThreeVector halfCircle6WrtBRCorner = shl6WrtBRCorner + G4ThreeVector(0.5*dp_shl6EmptyDimX,dp_resonatorAssemblyCurveCentralRadius,0.0);
  G4VPhysicalVolume * halfCircle6Empty = new G4PVPlacement(0,halfCircle6WrtBRCorner+brCornerOfBaseNbLayer,log_halfCircle6Empty,halfCircle6EmptyName,log_baseNbLayer,false,0,true);
  log_halfCircle6Empty->SetVisAttributes(air_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Vacuum",halfCircle6EmptyName,halfCircle6Empty));

  //------------------------------------------------------
  //HalfCircle 5 (conductor)
  G4String halfCircle6ConductorName = pName + "_halfCircle6Conductor";
  G4String halfCircle6ConductorNameSolid = halfCircle6ConductorName + "_solid";
  G4String halfCircle6ConductorNameLog = halfCircle6ConductorName + "_log";
  G4Tubs * solid_halfCircle6Conductor = new G4Tubs(halfCircle6ConductorNameSolid,dp_resonatorAssemblyCurveCentralRadius - dp_tlCouplingConductorDimY/2.0,dp_resonatorAssemblyCurveCentralRadius + dp_tlCouplingConductorDimY/2.0,dp_curveEmptyDimZ/2.0,270.*deg,180.*deg);
  G4LogicalVolume * log_halfCircle6Conductor = new G4LogicalVolume(solid_halfCircle6Conductor,niobium_mat,halfCircle6ConductorNameLog);  
  G4VPhysicalVolume * halfCircle6Conductor = new G4PVPlacement(0,G4ThreeVector(0,0,0),log_halfCircle6Conductor,halfCircle6ConductorName,log_halfCircle6Empty,false,0,true);
  log_halfCircle6Conductor->SetVisAttributes(niobium_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Niobium",halfCircle6ConductorName,halfCircle6Conductor));

  //------------------------------------------------------
  //Straight horizontal line (SHL) 7, (empty/cavity)
  G4String shl7EmptyName = pName + "_shl7Empty";
  G4String shl7EmptyNameSolid = shl7EmptyName + "_solid";
  G4String shl7EmptyNameLog = shl7EmptyName + "_log";
  G4Box * solid_shl7Empty = new G4Box(shl7EmptyNameSolid,0.5 * dp_shl7EmptyDimX,0.5 * dp_shl7EmptyDimY,0.5 * dp_shl7EmptyDimZ);
  G4LogicalVolume * log_shl7Empty = new G4LogicalVolume(solid_shl7Empty,air_mat,shl7EmptyNameLog);
  G4ThreeVector shl7WrtBRCorner = halfCircle6WrtBRCorner + G4ThreeVector(-0.5*dp_shl7EmptyDimX,dp_resonatorAssemblyCurveCentralRadius,0);
  G4VPhysicalVolume * shl7Empty = new G4PVPlacement(0,shl7WrtBRCorner+brCornerOfBaseNbLayer,log_shl7Empty,shl7EmptyName,log_baseNbLayer,false,0,true);
  log_shl7Empty->SetVisAttributes(air_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Vacuum",shl7EmptyName,shl7Empty));

  //------------------------------------------------------
  //Straight horizontal line (SHL) 7 (conductor)
  G4String shl7ConductorName = pName + "_shl7Conductor";
  G4String shl7ConductorNameSolid = shl7ConductorName + "_solid";
  G4String shl7ConductorNameLog = shl7ConductorName + "_log";
  G4Box * solid_shl7Conductor = new G4Box(shl7ConductorNameSolid,0.5 * dp_shl7ConductorDimX,0.5 * dp_shl7ConductorDimY,0.5 * dp_shl7ConductorDimZ);
  G4LogicalVolume * log_shl7Conductor = new G4LogicalVolume(solid_shl7Conductor,niobium_mat,shl7ConductorNameLog);
  G4VPhysicalVolume * shl7Conductor = new G4PVPlacement(0,G4ThreeVector(0,0,0),log_shl7Conductor,shl7ConductorName,log_shl7Empty,false,0,true);
  log_shl7Conductor->SetVisAttributes(niobium_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Niobium",shl7ConductorName,shl7Conductor));



  //------------------------------------------------------
  //Curve 3 (empty/cavity)
  G4String curve3EmptyName = pName + "_curve3Empty";
  G4String curve3EmptyNameSolid = curve3EmptyName + "_solid";
  G4String curve3EmptyNameLog = curve3EmptyName + "_log";
  G4Tubs * solid_curve3Empty = new G4Tubs(curve3EmptyNameSolid,dp_resonatorAssemblyCurveSmallestRadius,dp_resonatorAssemblyCurveSmallestRadius + dp_tlCouplingEmptyDimY,dp_curveEmptyDimZ/2.0,180.*deg,90.*deg);
  G4LogicalVolume * log_curve3Empty = new G4LogicalVolume(solid_curve3Empty,air_mat,curve3EmptyNameLog);  
  G4ThreeVector curve3WrtBRCorner = shl7WrtBRCorner + G4ThreeVector(-0.5*dp_shl7EmptyDimX,dp_resonatorAssemblyCurveCentralRadius,0.0); //Good for empty or conductor
  G4VPhysicalVolume * curve3Empty = new G4PVPlacement(0,curve3WrtBRCorner+brCornerOfBaseNbLayer,log_curve3Empty,curve3EmptyName,log_baseNbLayer,false,0,true);
  log_curve3Empty->SetVisAttributes(air_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Vacuum",curve3EmptyName,curve3Empty));


  //------------------------------------------------------
  //Curve 3 (conductor)
  G4String curve3ConductorName = pName + "_curve3Conductor";
  G4String curve3ConductorNameSolid = curve3ConductorName + "_solid";
  G4String curve3ConductorNameLog = curve3ConductorName + "_log";
  G4Tubs * solid_curve3Conductor = new G4Tubs(curve3ConductorNameSolid,dp_resonatorAssemblyCurveCentralRadius - dp_tlCouplingConductorDimY/2.0,dp_resonatorAssemblyCurveCentralRadius + dp_tlCouplingConductorDimY/2.0,dp_curveEmptyDimZ/2.0,180.*deg,90.*deg);
  G4LogicalVolume * log_curve3Conductor = new G4LogicalVolume(solid_curve3Conductor,niobium_mat,curve3ConductorNameLog);  
  G4VPhysicalVolume * curve3Conductor = new G4PVPlacement(0,G4ThreeVector(0,0,0),log_curve3Conductor,curve3ConductorName,log_curve3Empty,false,0,true);
  log_curve3Conductor->SetVisAttributes(niobium_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Niobium",curve3ConductorName,curve3Conductor));




  
  //------------------------------------------------------
  //Straight vertical line (SVL) 1, (empty/cavity)
  G4String svl1EmptyName = pName + "_svl1Empty";
  G4String svl1EmptyNameSolid = svl1EmptyName + "_solid";
  G4String svl1EmptyNameLog = svl1EmptyName + "_log";
  G4Box * solid_svl1Empty = new G4Box(svl1EmptyNameSolid,0.5 * dp_svl1EmptyDimX,0.5 * dp_svl1EmptyDimY,0.5 * dp_svl1EmptyDimZ);
  G4LogicalVolume * log_svl1Empty = new G4LogicalVolume(solid_svl1Empty,air_mat,svl1EmptyNameLog);
  G4ThreeVector svl1WrtBRCorner = curve3WrtBRCorner + G4ThreeVector(-1*dp_resonatorAssemblyCurveCentralRadius,0.5*dp_svl1EmptyDimY,0);
  G4VPhysicalVolume * svl1Empty = new G4PVPlacement(0,svl1WrtBRCorner+brCornerOfBaseNbLayer,log_svl1Empty,svl1EmptyName,log_baseNbLayer,false,0,true);
  log_svl1Empty->SetVisAttributes(air_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Vacuum",svl1EmptyName,svl1Empty));


  //------------------------------------------------------
  //Straight horizontal line (SHL) 7 (conductor)
  G4String svl1ConductorName = pName + "_svl1Conductor";
  G4String svl1ConductorNameSolid = svl1ConductorName + "_solid";
  G4String svl1ConductorNameLog = svl1ConductorName + "_log";
  G4Box * solid_svl1Conductor = new G4Box(svl1ConductorNameSolid,0.5 * dp_svl1ConductorDimX,0.5 * dp_svl1ConductorDimY,0.5 * dp_svl1ConductorDimZ);
  G4LogicalVolume * log_svl1Conductor = new G4LogicalVolume(solid_svl1Conductor,niobium_mat,svl1ConductorNameLog);
  G4VPhysicalVolume * svl1Conductor = new G4PVPlacement(0,G4ThreeVector(0,0,0),log_svl1Conductor,svl1ConductorName,log_svl1Empty,false,0,true);
  log_svl1Conductor->SetVisAttributes(niobium_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Niobium",svl1ConductorName,svl1Conductor));





  //------------------------------------------------------
  //Shunt coupler horizontal (empty/cavity)
  G4String shuntCouplerEmptyName = pName + "_shuntCouplerEmpty";
  G4String shuntCouplerEmptyNameSolid = shuntCouplerEmptyName + "_solid";
  G4String shuntCouplerEmptyNameLog = shuntCouplerEmptyName + "_log";
  G4Box * solid_shuntCouplerHorizontalEmpty = new G4Box("shuntCouplerHorizontalEmptySolid",0.5 * dp_shuntCouplerHorizontalEmptyDimX,0.5 * dp_shuntCouplerHorizontalEmptyDimY,0.5 * dp_shuntCouplerHorizontalEmptyDimZ);
  G4Box * solid_shuntCouplerLeftLobeEmpty = new G4Box("shuntCouplerLeftLobeEmpty",0.5 * dp_shuntCouplerLobeEmptyDimX,0.5* dp_shuntCouplerLobeEmptyDimY, 0.5* dp_shuntCouplerLobeEmptyDimZ);
  G4Box * solid_shuntCouplerRightLobeEmpty = new G4Box("shuntCouplerRightLobeEmpty",0.5 * dp_shuntCouplerLobeEmptyDimX,0.5* dp_shuntCouplerLobeEmptyDimY, 0.5* dp_shuntCouplerLobeEmptyDimZ);
  G4UnionSolid * solid_shuntCouplerMerge1Empty = new G4UnionSolid("shuntCouplerHorizontalPlusLeft",solid_shuntCouplerHorizontalEmpty,solid_shuntCouplerLeftLobeEmpty,0,G4ThreeVector(-0.5*dp_shuntCouplerHorizontalEmptyDimX + 0.5*dp_shuntCouplerLobeEmptyDimX,0.5*(dp_shuntCouplerHorizontalEmptyDimY+dp_shuntCouplerLobeEmptyDimY),0));
  G4UnionSolid * solid_shuntCouplerEmpty = new G4UnionSolid(shuntCouplerEmptyNameSolid,solid_shuntCouplerMerge1Empty,solid_shuntCouplerRightLobeEmpty,0,G4ThreeVector(0.5*dp_shuntCouplerHorizontalEmptyDimX - 0.5*dp_shuntCouplerLobeEmptyDimX,0.5*(dp_shuntCouplerHorizontalEmptyDimY+dp_shuntCouplerLobeEmptyDimY),0));    
  G4LogicalVolume * log_shuntCouplerEmpty = new G4LogicalVolume(solid_shuntCouplerEmpty,air_mat,shuntCouplerEmptyNameLog);  
  G4ThreeVector shuntCouplerEmptyWrtBRCorner = svl1WrtBRCorner + G4ThreeVector(0,0.5*(dp_shuntCouplerHorizontalEmptyDimY+dp_svl1EmptyDimY),0);
  G4VPhysicalVolume * shuntCouplerEmpty = new G4PVPlacement(0,shuntCouplerEmptyWrtBRCorner+brCornerOfBaseNbLayer,log_shuntCouplerEmpty,shuntCouplerEmptyName,log_baseNbLayer,false,0,true);
  log_shuntCouplerEmpty->SetVisAttributes(air_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Vacuum",shuntCouplerEmptyName,shuntCouplerEmpty));
  
  
  //------------------------------------------------------
  //Shunt coupler horizontal (conductor)
  G4String shuntCouplerConductorName = pName + "_shuntCouplerConductor";
  G4String shuntCouplerConductorNameSolid = shuntCouplerConductorName + "_solid";
  G4String shuntCouplerConductorNameLog = shuntCouplerConductorName + "_log";
  G4Box * solid_shuntCouplerHorizontalConductorSansNub = new G4Box("ShuntCouplerHorizontalConductorSansNub",0.5 * dp_shuntCouplerHorizontalConductorDimX,0.5 * dp_shuntCouplerHorizontalConductorDimY,0.5 * dp_shuntCouplerHorizontalConductorDimZ);
  G4Box * solid_shuntCouplerHorizontalConductorNub = new G4Box("ShuntCouplerHorizontalConductorNub",0.5 * dp_shuntCouplerHorizontalConductorNubDimX,0.5 * dp_shuntCouplerHorizontalConductorNubDimY,0.5 * dp_shuntCouplerHorizontalConductorNubDimZ);
  G4UnionSolid * solid_shuntCouplerHorizontalConductor = new G4UnionSolid("shuntCouplerHorizontalConductor",solid_shuntCouplerHorizontalConductorSansNub,solid_shuntCouplerHorizontalConductorNub,0,G4ThreeVector(0,-0.5*(dp_shuntCouplerHorizontalConductorDimY+dp_shuntCouplerHorizontalConductorNubDimY)));
  G4Box * solid_shuntCouplerLeftLobeConductor = new G4Box("shuntCouplerLeftLobeConductor",0.5 * dp_shuntCouplerLobeConductorDimX,0.5* dp_shuntCouplerLobeConductorDimY, 0.5* dp_shuntCouplerLobeConductorDimZ);
  G4Box * solid_shuntCouplerRightLobeConductor = new G4Box("shuntCouplerRightLobeConductor",0.5 * dp_shuntCouplerLobeConductorDimX,0.5* dp_shuntCouplerLobeConductorDimY, 0.5* dp_shuntCouplerLobeConductorDimZ);
  G4UnionSolid * solid_shuntCouplerMerge2Conductor = new G4UnionSolid("shuntCouplerHorizontalPlusNubPlusLeft",solid_shuntCouplerHorizontalConductor,solid_shuntCouplerLeftLobeConductor,0,G4ThreeVector(-0.5*dp_shuntCouplerHorizontalConductorDimX + 0.5*dp_shuntCouplerLobeConductorDimX,0.5*(dp_shuntCouplerHorizontalConductorDimY+dp_shuntCouplerLobeConductorDimY),0));
  G4UnionSolid * solid_shuntCouplerConductor = new G4UnionSolid(shuntCouplerConductorNameSolid,solid_shuntCouplerMerge2Conductor,solid_shuntCouplerRightLobeConductor,0,G4ThreeVector(0.5*dp_shuntCouplerHorizontalConductorDimX - 0.5*dp_shuntCouplerLobeConductorDimX,0.5*(dp_shuntCouplerHorizontalConductorDimY+dp_shuntCouplerLobeConductorDimY),0));      
  G4LogicalVolume * log_shuntCouplerConductor = new G4LogicalVolume(solid_shuntCouplerConductor,niobium_mat,shuntCouplerConductorNameLog);
  G4VPhysicalVolume * shuntCouplerConductor = new G4PVPlacement(0,G4ThreeVector(0,0,0),log_shuntCouplerConductor,shuntCouplerConductorName,log_shuntCouplerEmpty,false,0,true);
  log_shuntCouplerConductor->SetVisAttributes(niobium_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Niobium",shuntCouplerConductorName,shuntCouplerConductor));
 
  
									  
  
  
  
}
  


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<std::tuple<std::string,G4String,G4VPhysicalVolume*> > RISQTutorialResonatorAssembly::GetListOfAllFundamentalSubVolumes()
{
  return fFundamentalVolumeList;
}

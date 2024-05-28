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
#include "RISQTutorialCornerFluxLine.hh"
#include "RISQTutorialPad.hh"
#include "RISQTutorialDetectorParameters.hh"

using namespace RISQTutorialDetectorParameters;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Primary Constructor
RISQTutorialCornerFluxLine::RISQTutorialCornerFluxLine(G4RotationMatrix * pRot,
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

  ConstructCornerFluxLine(pRot,
			    tLate,
			    pName,
			    pMotherLogical,
			    pMany,
			    pCopyNo,
			    pSurfChk);

  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Default Constructor
RISQTutorialCornerFluxLine::RISQTutorialCornerFluxLine()
{  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Destructor
RISQTutorialCornerFluxLine::~RISQTutorialCornerFluxLine()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//Moving implementation down here so it's not in the constructor
void RISQTutorialCornerFluxLine::ConstructCornerFluxLine(G4RotationMatrix * pRot,
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

  //Here, we can just make our baseNbLayer a simple box, since it won't interfere/overlap with any other
  //elements at this layer in the heirarchy
  G4Box * solid_baseNbLayer = new G4Box(baseNbLayerNameSolid,0.5*dp_cornerFluxLineBaseNbLayerDimX,0.5*dp_cornerFluxLineBaseNbLayerDimY,0.5*dp_cornerFluxLineBaseNbLayerDimZ);

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
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Niobium",baseNbLayerName,phys_baseNbLayer));



  //------------------------------------------------------------------------------------------
  //Now make the pad of the flux line. The logical mother volume of these is the base
  //Nb layer. Pads come with their own visualization attributes already set.
  
  //Pad 1:
  double cornerFluxLinePadCenterOffsetFromTopOrSide = pow(pow((0.5*dp_padEmptyPart1DimX),2) + pow(0.5*dp_padEmptyPart1DimY,2),0.5);//Pow is not const...
  double cornerFluxLinePadOffsetY = 0.5*dp_cornerFluxLineBaseNbLayerDimY - cornerFluxLinePadCenterOffsetFromTopOrSide;
  double cornerFluxLinePadOffsetX = -0.5*dp_cornerFluxLineBaseNbLayerDimX + cornerFluxLinePadCenterOffsetFromTopOrSide;
  G4String pad1Name = pName + "_FluxLinePad1";
  G4RotationMatrix * pad1Rot = new G4RotationMatrix();
  pad1Rot->rotateZ(45.*deg);
  RISQTutorialPad * pad1 = new RISQTutorialPad(pad1Rot,
							       G4ThreeVector(cornerFluxLinePadOffsetX,cornerFluxLinePadOffsetY,0),
							       pad1Name,
							       log_baseNbLayer,
							       false,
							       0,
							       checkOverlaps);
  G4LogicalVolume * log_pad1 = pad1->GetLogicalVolume();
  G4VPhysicalVolume * phys_pad1 = pad1->GetPhysicalVolume();
  AddComplexGeometryPadSubVolumesToThisList(pad1);
  

  //------------------------------------------------------------------------------------------
  //Now make the flux line itself, in several parts:
  //1. Curved part emerging from pad
  //2. Straight horizontal
  //3. Curved part down
  //4. Straight vertical
  //Each of these has empty and conductor elements

  //Some mathematical things to help (can't put these in the dp file because they can't be executed as constexpr)
  double cornerFluxLineCurve1WrtBaseLayerX = dp_cornerFluxLineCurve1WrtPadCenterX + cornerFluxLinePadOffsetX;
  double cornerFluxLineCurve1WrtBaseLayerY = dp_cornerFluxLineCurve1WrtPadCenterY + cornerFluxLinePadOffsetY;

  
  //Curve 1, empty part
  G4String curve1NameEmpty = pName + "_curve1Empty";
  G4String curve1NameEmptySolid = curve1NameEmpty + "_solid";
  G4String curve1NameEmptyLog = curve1NameEmpty + "_log";
  G4Tubs * solid_curve1Empty = new G4Tubs(curve1NameEmptySolid,dp_cornerFluxLineCurveRadius - dp_fluxLineEmptyDimX/2.0, dp_cornerFluxLineCurveRadius + dp_fluxLineEmptyDimX/2.0,dp_curveEmptyDimZ/2.0,225.0*deg,45.0*deg);

  G4LogicalVolume * log_curve1Empty = new G4LogicalVolume(solid_curve1Empty,air_mat,curve1NameEmptyLog);  
  G4ThreeVector curve1EmptyWrtBaseLayerCenter(cornerFluxLineCurve1WrtBaseLayerX,cornerFluxLineCurve1WrtBaseLayerY,0.0);
  G4VPhysicalVolume * curve1Empty = new G4PVPlacement(0,curve1EmptyWrtBaseLayerCenter,log_curve1Empty,curve1NameEmpty,log_baseNbLayer,false,0,true);
  log_curve1Empty->SetVisAttributes(air_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Vacuum",curve1NameEmpty,curve1Empty));


  //Curve 1, conductor part
  G4String curve1NameConductor = pName + "_curve1Conductor";
  G4String curve1NameConductorSolid = curve1NameConductor + "_solid";
  G4String curve1NameConductorLog = curve1NameConductor + "_log";
  G4Tubs * solid_curve1Conductor = new G4Tubs(curve1NameConductorSolid,dp_cornerFluxLineCurveRadius - dp_fluxLineConductorDimX/2.0, dp_cornerFluxLineCurveRadius + dp_fluxLineConductorDimX/2.0,dp_curveEmptyDimZ/2.0,225.0*deg,45.0*deg);

  G4LogicalVolume * log_curve1Conductor = new G4LogicalVolume(solid_curve1Conductor,niobium_mat,curve1NameConductorLog);  
  G4VPhysicalVolume * curve1Conductor = new G4PVPlacement(0,G4ThreeVector(0,0,0),log_curve1Conductor,curve1NameConductor,log_curve1Empty,false,0,true);
  log_curve1Conductor->SetVisAttributes(niobium_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Niobium",curve1NameConductor,curve1Conductor));



  //Horizontal line, empty part
  G4String horizontalNameEmpty = pName + "_horizontalEmpty";
  G4String horizontalNameEmptySolid = horizontalNameEmpty + "_solid";
  G4String horizontalNameEmptyLog = horizontalNameEmpty + "_log";
  G4Box * solid_horizontalEmpty = new G4Box(horizontalNameEmptySolid,0.5*dp_cornerFluxLineHorizontalEmptyDimX,0.5*dp_cornerFluxLineHorizontalEmptyDimY,0.5*dp_cornerFluxLineHorizontalEmptyDimZ);

  G4LogicalVolume * log_horizontalEmpty = new G4LogicalVolume(solid_horizontalEmpty,air_mat,horizontalNameEmptyLog);  
  G4ThreeVector horizontalEmptyWrtCurve1Empty(0.5*dp_cornerFluxLineHorizontalEmptyDimX,-dp_cornerFluxLineCurveRadius,0.0);
  G4VPhysicalVolume * horizontalEmpty = new G4PVPlacement(0,horizontalEmptyWrtCurve1Empty+curve1EmptyWrtBaseLayerCenter,log_horizontalEmpty,horizontalNameEmpty,log_baseNbLayer,false,0,true);
  log_horizontalEmpty->SetVisAttributes(air_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Vacuum",horizontalNameEmpty,horizontalEmpty));


  //Horizontal line, conductor part
  G4String horizontalNameConductor = pName + "_horizontalConductor";
  G4String horizontalNameConductorSolid = horizontalNameConductor + "_solid";
  G4String horizontalNameConductorLog = horizontalNameConductor + "_log";
  G4Box * solid_horizontalConductor = new G4Box(horizontalNameConductorSolid,0.5*dp_cornerFluxLineHorizontalConductorDimX,0.5*dp_cornerFluxLineHorizontalConductorDimY,0.5*dp_cornerFluxLineHorizontalConductorDimZ);

  G4LogicalVolume * log_horizontalConductor = new G4LogicalVolume(solid_horizontalConductor,niobium_mat,horizontalNameConductorLog);  
  G4VPhysicalVolume * horizontalConductor = new G4PVPlacement(0,G4ThreeVector(0,0,0),log_horizontalConductor,horizontalNameConductor,log_horizontalEmpty,false,0,true);
  log_horizontalConductor->SetVisAttributes(niobium_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Niobium",horizontalNameConductor,horizontalConductor));


  //Curve 2, empty part
  G4String curve2NameEmpty = pName + "_curve2Empty";
  G4String curve2NameEmptySolid = curve2NameEmpty + "_solid";
  G4String curve2NameEmptyLog = curve2NameEmpty + "_log";
  G4Tubs * solid_curve2Empty = new G4Tubs(curve2NameEmptySolid,dp_cornerFluxLineCurveRadius - dp_fluxLineEmptyDimX/2.0, dp_cornerFluxLineCurveRadius + dp_fluxLineEmptyDimX/2.0,dp_curveEmptyDimZ/2.0,0.0*deg,90.0*deg);

  G4LogicalVolume * log_curve2Empty = new G4LogicalVolume(solid_curve2Empty,air_mat,curve2NameEmptyLog);  
  G4ThreeVector curve2EmptyWrtHorizontalEmpty(0.5*dp_cornerFluxLineHorizontalEmptyDimX,-dp_cornerFluxLineCurveRadius,0.0);
  G4VPhysicalVolume * curve2Empty = new G4PVPlacement(0,curve2EmptyWrtHorizontalEmpty+horizontalEmptyWrtCurve1Empty+curve1EmptyWrtBaseLayerCenter,log_curve2Empty,curve2NameEmpty,log_baseNbLayer,false,0,true);
  log_curve2Empty->SetVisAttributes(air_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Vacuum",curve2NameEmpty,curve2Empty));



  //Curve 2, conductor part
  G4String curve2NameConductor = pName + "_curve2Conductor";
  G4String curve2NameConductorSolid = curve2NameConductor + "_solid";
  G4String curve2NameConductorLog = curve2NameConductor + "_log";
  G4Tubs * solid_curve2Conductor = new G4Tubs(curve2NameConductorSolid,dp_cornerFluxLineCurveRadius - dp_fluxLineConductorDimX/2.0, dp_cornerFluxLineCurveRadius + dp_fluxLineConductorDimX/2.0,dp_curveEmptyDimZ/2.0,0.0*deg,90.0*deg);

  G4LogicalVolume * log_curve2Conductor = new G4LogicalVolume(solid_curve2Conductor,niobium_mat,curve2NameConductorLog);  
  G4VPhysicalVolume * curve2Conductor = new G4PVPlacement(0,G4ThreeVector(0,0,0),log_curve2Conductor,curve2NameConductor,log_curve2Empty,false,0,true);
  log_curve2Conductor->SetVisAttributes(niobium_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Niobium",curve2NameConductor,curve2Conductor));





  //Vertical line, empty part
  G4String verticalNameEmpty = pName + "_verticalEmpty";
  G4String verticalNameEmptySolid = verticalNameEmpty + "_solid";
  G4String verticalNameEmptyLog = verticalNameEmpty + "_log";
  G4Box * solid_verticalEmpty = new G4Box(verticalNameEmptySolid,0.5*dp_cornerFluxLineVerticalEmptyDimX,0.5*dp_cornerFluxLineVerticalEmptyDimY,0.5*dp_cornerFluxLineVerticalEmptyDimZ);

  G4LogicalVolume * log_verticalEmpty = new G4LogicalVolume(solid_verticalEmpty,air_mat,verticalNameEmptyLog);  
  G4ThreeVector verticalEmptyWrtCurve2Empty(dp_cornerFluxLineCurveRadius,-0.5*dp_cornerFluxLineVerticalEmptyDimY,0.0);
  G4VPhysicalVolume * verticalEmpty = new G4PVPlacement(0,verticalEmptyWrtCurve2Empty+curve2EmptyWrtHorizontalEmpty+horizontalEmptyWrtCurve1Empty+curve1EmptyWrtBaseLayerCenter,log_verticalEmpty,verticalNameEmpty,log_baseNbLayer,false,0,true);
  log_verticalEmpty->SetVisAttributes(air_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Vacuum",verticalNameEmpty,verticalEmpty));

  /*
  std::cout << "---...............^^^^^" << std::endl;
  verticalEmpty->CheckOverlaps(1000000,0,true);
  std::cout << "---...............^^^^^" << std::endl;
  */
  
  

  //Vertical line, conductor part
  G4String verticalNameConductor = pName + "_verticalConductor";
  G4String verticalNameConductorSolid = verticalNameConductor + "_solid";
  G4String verticalNameConductorLog = verticalNameConductor + "_log";
  G4Box * solid_verticalConductor = new G4Box(verticalNameConductorSolid,0.5*dp_cornerFluxLineVerticalConductorDimX,0.5*dp_cornerFluxLineVerticalConductorDimY,0.5*dp_cornerFluxLineVerticalConductorDimZ);

  G4LogicalVolume * log_verticalConductor = new G4LogicalVolume(solid_verticalConductor,niobium_mat,verticalNameConductorLog);  
  G4VPhysicalVolume * verticalConductor = new G4PVPlacement(0,G4ThreeVector(0,0,0),log_verticalConductor,verticalNameConductor,log_verticalEmpty,false,0,true);
  log_verticalConductor->SetVisAttributes(niobium_vis);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Niobium",verticalNameConductor,verticalConductor));













  


  /*

  G4String curve1EmptyName = pName + "_curve1Empty";
  G4String curve1EmptyNameSolid = curve1EmptyName + "_solid";
  G4String curve1EmptyNameLog = curve1EmptyName + "_log";
  G4Tubs * solid_curve1Empty = new G4Tubs(curve1EmptyNameSolid,dp_resonatorAssemblyCurveSmallestRadius,dp_resonatorAssemblyCurveSmallestRadius + dp_tlCouplingEmptyDimY,dp_curveEmptyDimZ/2.0,180.*deg,90.*deg);
  G4LogicalVolume * log_curve1Empty = new G4LogicalVolume(solid_curve1Empty,air_mat,curve1EmptyNameLog);  
  G4ThreeVector curve1WrtBRCorner(-1*dp_tlCouplingEmptyDimX,0.5*dp_tlCouplingEmptyDimY + dp_resonatorAssemblyBaseNbEdgeBottomDimY + dp_resonatorAssemblyCurveCentralRadius,0.0); //Good for empty or conductor
  G4VPhysicalVolume * curve1Empty = new G4PVPlacement(0,curve1WrtBRCorner+brCornerOfBaseNbLayer,log_curve1Empty,curve1EmptyName,log_baseNbLayer,false,0,true);
  log_curve1Empty->SetVisAttributes(air_vis);



  G4Box * solid_fluxLineEmpty = new G4Box(flNameEmptySolid,
					  0.5 * dp_fluxLineEmptyDimX,
					  0.5 * dp_fluxLineEmptyDimY,
					  0.5 * dp_fluxLineEmptyDimZ);
  
  G4LogicalVolume * log_fluxLineEmpty = new G4LogicalVolume(solid_fluxLineEmpty,
							    air_mat,
							    flNameEmptyLog);
  G4VPhysicalVolume * phys_fluxLineEmpty = new G4PVPlacement(0,
							     G4ThreeVector(0,dp_fluxLineLineOffsetY,0),
							     log_fluxLineEmpty,
							     flNameEmpty,
							     log_baseNbLayer,
							     false,
							     0,
							     true);
  log_fluxLineEmpty->SetVisAttributes(air_vis);
  
  
  G4String flNameConductor = pName + "_FluxLineConductor";
  G4String flNameConductorSolid = flNameConductor + "_solid";
  G4String flNameConductorLog = flNameConductor + "_log";
  G4Box * solid_fluxLineConductor = new G4Box(flNameConductorSolid,
					      0.5 * dp_fluxLineConductorDimX,
					      0.5 * dp_fluxLineConductorDimY,
					      0.5 * dp_fluxLineConductorDimZ);  
  G4LogicalVolume * log_fluxLineConductor = new G4LogicalVolume(solid_fluxLineConductor,
								niobium_mat,
								flNameConductorLog);
  G4VPhysicalVolume * phys_fluxLineConductor = new G4PVPlacement(0,
								 G4ThreeVector(0,0,0),
								 log_fluxLineConductor,
								 flNameConductor,
								 log_fluxLineEmpty,
								 false,
								 0,
								 true);
  log_fluxLineConductor->SetVisAttributes(niobium_vis);  
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
std::vector<std::tuple<std::string,G4String,G4VPhysicalVolume*> > RISQTutorialCornerFluxLine::GetListOfAllFundamentalSubVolumes()
{
  return fFundamentalVolumeList;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RISQTutorialCornerFluxLine::AddComplexGeometryPadSubVolumesToThisList(RISQTutorialPad * pad)
{
  for( int iSubVol = 0; iSubVol < pad->GetListOfAllFundamentalSubVolumes().size(); ++iSubVol){
    std::tuple<std::string,G4String,G4VPhysicalVolume*> theTuple(std::get<0>(pad->GetListOfAllFundamentalSubVolumes()[iSubVol]),
								 std::get<1>(pad->GetListOfAllFundamentalSubVolumes()[iSubVol]),
								 std::get<2>(pad->GetListOfAllFundamentalSubVolumes()[iSubVol]));
    fFundamentalVolumeList.push_back(theTuple);
  }
}


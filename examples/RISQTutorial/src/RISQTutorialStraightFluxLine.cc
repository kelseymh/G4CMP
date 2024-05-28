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
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

//Includes (specific to this project)
#include "RISQTutorialStraightFluxLine.hh"
#include "RISQTutorialPad.hh"
#include "RISQTutorialDetectorParameters.hh"

using namespace RISQTutorialDetectorParameters;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Primary Constructor
RISQTutorialStraightFluxLine::RISQTutorialStraightFluxLine(G4RotationMatrix * pRot,
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

  ConstructStraightFluxLine(pRot,
			    tLate,
			    pName,
			    pMotherLogical,
			    pMany,
			    pCopyNo,
			    pSurfChk);

  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Default Constructor
RISQTutorialStraightFluxLine::RISQTutorialStraightFluxLine()
{  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Destructor
RISQTutorialStraightFluxLine::~RISQTutorialStraightFluxLine()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//Moving implementation down here so it's not in the constructor
void RISQTutorialStraightFluxLine::ConstructStraightFluxLine(G4RotationMatrix * pRot,
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
  G4Box * solid_baseNbLayer = new G4Box(baseNbLayerNameSolid,0.5*dp_fluxLineBaseNbLayerDimX,0.5*dp_fluxLineBaseNbLayerDimY,0.5*dp_fluxLineBaseNbLayerDimZ);

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
  //Ni layer. Pads come with their own visualization attributes already set.
  
  //Pad 1:
  G4String pad1Name = pName + "_FluxLinePad1";
  G4RotationMatrix * pad1Rot = new G4RotationMatrix();
  pad1Rot->rotateZ(90.*deg);
  RISQTutorialPad * pad1 = new RISQTutorialPad(pad1Rot,
							       G4ThreeVector(0,dp_fluxLinePadOffsetY,0),
							       pad1Name,
							       log_baseNbLayer,
							       false,
							       0,
							       checkOverlaps);
  G4LogicalVolume * log_pad1 = pad1->GetLogicalVolume();
  G4VPhysicalVolume * phys_pad1 = pad1->GetPhysicalVolume();
  AddComplexGeometryPadSubVolumesToThisList(pad1);
  
  

  //------------------------------------------------------------------------------------------
  //Now make the flux line itself, in two parts: empty and conductor. The logical
  //mother volume of this is the base Ni layer.
  G4String flNameEmpty = pName + "_FluxLineEmpty";
  G4String flNameEmptySolid = flNameEmpty + "_solid";
  G4String flNameEmptyLog = flNameEmpty + "_log";
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
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Vacuum",flNameEmpty,phys_fluxLineEmpty));
  
  
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
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Niobium",flNameConductor,phys_fluxLineConductor));
  
  



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
std::vector<std::tuple<std::string,G4String,G4VPhysicalVolume*> > RISQTutorialStraightFluxLine::GetListOfAllFundamentalSubVolumes()
{
  return fFundamentalVolumeList;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RISQTutorialStraightFluxLine::AddComplexGeometryPadSubVolumesToThisList(RISQTutorialPad * pad)
{
  for( int iSubVol = 0; iSubVol < pad->GetListOfAllFundamentalSubVolumes().size(); ++iSubVol){
    std::tuple<std::string,G4String,G4VPhysicalVolume*> theTuple(std::get<0>(pad->GetListOfAllFundamentalSubVolumes()[iSubVol]),
								 std::get<1>(pad->GetListOfAllFundamentalSubVolumes()[iSubVol]),
								 std::get<2>(pad->GetListOfAllFundamentalSubVolumes()[iSubVol]));
    fFundamentalVolumeList.push_back(theTuple);
  }
}


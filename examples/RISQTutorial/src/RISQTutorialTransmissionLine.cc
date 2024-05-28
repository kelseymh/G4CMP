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
#include "RISQTutorialTransmissionLine.hh"
#include "RISQTutorialPad.hh"
#include "RISQTutorialDetectorParameters.hh"

using namespace RISQTutorialDetectorParameters;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Primary Constructor
RISQTutorialTransmissionLine::RISQTutorialTransmissionLine(G4RotationMatrix * pRot,
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

  ConstructTransmissionLine(pRot,
			    tLate,
			    pName,
			    pMotherLogical,
			    pMany,
			    pCopyNo,
			    pSurfChk);

  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Default Constructor
RISQTutorialTransmissionLine::RISQTutorialTransmissionLine()
{  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Destructor
RISQTutorialTransmissionLine::~RISQTutorialTransmissionLine()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//Moving implementation down here so it's not in the constructor
void RISQTutorialTransmissionLine::ConstructTransmissionLine(G4RotationMatrix * pRot,
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

  //Unfortunately, we need our baseNbLayer to not interfere with the resonator structures. So we're going
  //to have to form a g4union object from the external pad shapes and use this as our base layer.
  G4UnionSolid * solid_baseNbLayer = CreatePieceBasedNbLayer(baseNbLayerNameSolid);

  
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
  
  //Push this sub volume (the niobium base layer) back into the fundamental volume list
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Niobium",baseNbLayerName,phys_baseNbLayer));
  



  //------------------------------------------------------------------------------------------
  //Now make the pads of the transmission line. The logical mother volume of these is the base
  //Nb layer. Pads come with their own visualization attributes already set.
  
  //Pad 1:
  G4String pad1Name = pName + "_TransmissionLinePad1";
  RISQTutorialPad * pad1 = new RISQTutorialPad(0,
							       G4ThreeVector(dp_transmissionLinePad1Offset,0,0),
							       pad1Name,
							       log_baseNbLayer,
							       false,
							       0,
							       checkOverlaps);
  G4LogicalVolume * log_pad1 = pad1->GetLogicalVolume();
  G4VPhysicalVolume * phys_pad1 = pad1->GetPhysicalVolume();

  //Loop through the fundamental sub-volumes and push them back into the fundamental subvolume list for the transmission line. We have to do this here
  //because the pads are composite volumes
  AddComplexGeometryPadSubVolumesToThisList(pad1);
  

  //Pad 2: rotate around Z axis by 180 degrees
  G4String pad2Name = pName + "_TransmissionLinePad2";
  G4RotationMatrix * pad2Rot = new G4RotationMatrix();
  pad2Rot->rotateZ(180*deg);
  RISQTutorialPad * pad2 = new RISQTutorialPad(pad2Rot,
							       G4ThreeVector(dp_transmissionLinePad2Offset,0,0),
							       pad2Name,
							       log_baseNbLayer,
							       false,
							       0,
							       checkOverlaps);
  G4LogicalVolume * log_pad2 = pad1->GetLogicalVolume();
  G4VPhysicalVolume * phys_pad2 = pad1->GetPhysicalVolume();
  
  //Loop through the fundamental sub-volumes and push them back into the fundamental subvolume list for the transmission line. We have to do this here
  //because the pads are composite volumes
  AddComplexGeometryPadSubVolumesToThisList(pad2);
  
  


  
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
								     log_baseNbLayer,
								     false,
								     0,
								     true);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Vacuum",tlNameEmpty,phys_transmissionLineEmpty));



  
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
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Niobium",tlNameConductor,phys_transmissionLineConductor));

  



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
//Used for building the base layer for the Niobium. Can't be a pure rectangle because that rectangle
//will overlap with the resonators, which we can't have... Ripping some of this from the code in
//the pad file and the transmission line code, so be careful if you end up having to change the dimensions of things...
G4UnionSolid * RISQTutorialTransmissionLine::CreatePieceBasedNbLayer(G4String nameSolid)
{
  //------------------------------------------------------------------------------------------
  //For the pad, we start with a volume of air with the same dimensions as the niobium substrate.
  //Then we add the niobium pad within.
  G4Box * solid_padEmptyPart1 = new G4Box("baseNbLayerEmptyPadPart1Solid",0.5 * dp_padEmptyPart1DimX,0.5 * dp_padEmptyPart1DimY,0.5 * dp_padEmptyPart1DimZ);
  G4Trd * solid_padEmptyPart2 = new G4Trd("baseNbLayerEmptyPadPart2Solid",0.5 * dp_padEmptyPart2TrdX1,0.5 * dp_padEmptyPart2TrdX2,0.5 * dp_padEmptyPart2TrdY1,0.5 * dp_padEmptyPart2TrdY2,0.5 * dp_padEmptyPart2TrdZ);

  //We need to rotate the part 2 so that it can be aligned and placed next to part 1
  G4RotationMatrix * rotEmptyPart2 = new G4RotationMatrix();
  rotEmptyPart2->rotateX(90.*deg);
  rotEmptyPart2->rotateY(-90.*deg);

  G4ThreeVector transPart2EmptyWrtPart1Empty(dp_padEmptyPart1DimX/2.0 + dp_padEmptyPart2TrdZ/2.0,0,0);
  G4UnionSolid * solid_padEmpty = new G4UnionSolid("BaseNbLayerEmptyPadSolid",
						   solid_padEmptyPart1,
						   solid_padEmptyPart2,
						   rotEmptyPart2,
						   transPart2EmptyWrtPart1Empty);
  
  G4Box * solid_transmissionLineEmpty = new G4Box("BaseNbLayerTransmissionLineSolid",
						  0.5 * (dp_transmissionLinePad2Offset - dp_transmissionLinePad1Offset - 2 * dp_padEmptyPart2TrdZ - 2 * 0.5 * dp_padEmptyPart1DimX),
						  0.5 * dp_transmissionLineCavityFullWidth,
						  0.5 * dp_transmissionLineBaseLayerDimZ);

  G4UnionSolid * solid_merger1 = new G4UnionSolid("Merger1",
						  solid_transmissionLineEmpty,
						  solid_padEmpty,
						  0,
						  G4ThreeVector(dp_transmissionLinePad1Offset,0,0));

  G4RotationMatrix * pad2Rot = new G4RotationMatrix();
  pad2Rot->rotateZ(180*deg);
  G4UnionSolid * solid_baseNbLayer = new G4UnionSolid(nameSolid,
						      solid_merger1,
						      solid_padEmpty,
						      pad2Rot,
						      G4ThreeVector(dp_transmissionLinePad2Offset,0,0));

  return solid_baseNbLayer;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<std::tuple<std::string,G4String,G4VPhysicalVolume*> > RISQTutorialTransmissionLine::GetListOfAllFundamentalSubVolumes()
{
  return fFundamentalVolumeList;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RISQTutorialTransmissionLine::AddComplexGeometryPadSubVolumesToThisList(RISQTutorialPad * pad)
{
  for( int iSubVol = 0; iSubVol < pad->GetListOfAllFundamentalSubVolumes().size(); ++iSubVol){
    std::tuple<std::string,G4String,G4VPhysicalVolume*> theTuple(std::get<0>(pad->GetListOfAllFundamentalSubVolumes()[iSubVol]),
								 std::get<1>(pad->GetListOfAllFundamentalSubVolumes()[iSubVol]),
								 std::get<2>(pad->GetListOfAllFundamentalSubVolumes()[iSubVol]));
    fFundamentalVolumeList.push_back(theTuple);
  }
}

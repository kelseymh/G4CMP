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
// 20260109  M. Kelsey -- G4CMP-569: Remove unused local variables.

/// \file ValidationTransmissionLine.cc
/// \brief Class implementation for the transmission line in the validation
///		example.

//Includes (basic)
#include "ValidationTransmissionLine.hh"
#include "G4LatticePhysical.hh"
#include "G4CMPLogicalBorderSurface.hh"
#include "ValidationDetectorParameters.hh"
#include "ValidationPad.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4Orb.hh"
#include "G4PVPlacement.hh"
#include "G4RunManager.hh"
#include "G4Sphere.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4Trd.hh"
#include "G4UnionSolid.hh"
#include "G4VisAttributes.hh"

using namespace ValidationDetectorParameters;

// Primary Constructor
ValidationTransmissionLine::
ValidationTransmissionLine(G4RotationMatrix * pRot,
                           const G4ThreeVector & tLate,
                           const G4String & pName,
                           G4LogicalVolume * pMotherLogical,
                           G4bool pMany,
                           G4int pCopyNo,
                           G4LatticeManager * LM,
                           std::map<std::string,G4LatticeLogical*> logicalLatticeContainer,
                           std::map<std::string,G4CMPSurfaceProperty*> borderContainer,
                           G4bool pSurfChk)

{
  //Here, use the inputs to this to set up the geometry and fill out the
  //PVPlacement data member, which is the real output from this class (and
  //which we'll access in our detector construction file.)

  ConstructTransmissionLine(pRot,tLate,pName,pMotherLogical,pMany,pCopyNo,LM,
                            logicalLatticeContainer,borderContainer,pSurfChk);



  
}


// Default Constructor
ValidationTransmissionLine::ValidationTransmissionLine() {  
}

// Destructor
ValidationTransmissionLine::~ValidationTransmissionLine() {
}

//Moving implementation down here so it's not in the constructor
void ValidationTransmissionLine::
ConstructTransmissionLine(G4RotationMatrix * pRot,const G4ThreeVector & tLate,
                          const G4String & pName,
                          G4LogicalVolume * pMotherLogical,G4bool pMany,
                          G4int pCopyNo,G4LatticeManager * LM,
                          std::map<std::string,G4LatticeLogical*> logicalLatticeContainer,
                          std::map<std::string,G4CMPSurfaceProperty*> borderContainer,
                          G4bool pSurfChk) {
  
  //Start with some preliminaries - NIST manager
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* aluminum_mat = nist->FindOrBuildMaterial("G4_Al");
  G4Material* air_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  //Set up lattice information
  if (logicalLatticeContainer.count("Aluminum") == 0) {
    std::cout << "Uh oh! Trying to access logicalLatticeContainer[Aluminum] "
              << "but it's not there..." << std::endl;
  }
  G4LatticeLogical* AlLogical = logicalLatticeContainer["Aluminum"];
  
  //Set up the aluminum visualization
  G4VisAttributes* aluminum_vis= new G4VisAttributes(G4Colour(0.0,1.0,1.0,0.5));
  aluminum_vis->SetVisibility(true);
  G4VisAttributes* air_vis= new G4VisAttributes(G4Colour(0.5,0.5,0.5,0.5));
  air_vis->SetVisibility(true);
  

  
  
  //Start with a base layer of aluminum into which our objects will fit. We'll
  //return this in the end.
  G4String baseAlLayerName = pName;
  G4String baseAlLayerNameSolid = pName + "_solid";
  G4String baseAlLayerNameLog = pName + "_log";

  //Unfortunately, we need our baseAlLayer to not interfere with the resonator
  //structures. So we're going to have to form a g4union object from the
  //external pad shapes and use this as our base layer.
  G4UnionSolid * solid_baseAlLayer =
    CreatePieceBasedAlLayer(baseAlLayerNameSolid);

  
  //Now attribute a physical material to the housing
  G4LogicalVolume * log_baseAlLayer
    = new G4LogicalVolume(solid_baseAlLayer,aluminum_mat,baseAlLayerNameLog);
  log_baseAlLayer->SetVisAttributes(G4VisAttributes::Invisible);

  //Now, create a physical volume and G4PVPlacement for storing as the final
  //output. This is the top volume.
  G4VPhysicalVolume* phys_baseAlLayer =
    new G4PVPlacement(pRot,tLate,log_baseAlLayer,baseAlLayerName,pMotherLogical,
                      pMany,pCopyNo,pSurfChk);


  //This base layer must be given an external interface with the
  //layer into which it's embedded (pMotherLogical), which requires a physical
  //volume, not passed into this code. We therefore add this to the list of
  //the toBeBoundariedVolumeList with info on what kind of boundary it must
  //have. Since this base layer is actually *fully* occupied by empty pad/TL
  //layers that come later, we will make it the same boundary as those (?
  //What happens there?)

  //Push this sub volume (the aluminum base layer) back into the fundamental
  //volume list
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Aluminum",baseAlLayerName,phys_baseAlLayer));
  
  //Now make the pads of the transmission line. The logical mother volume of
  //these is the base Al layer. Pads come with their own visualization
  //attributes already set.
  
  //Pad 1:
  G4String pad1Name = pName + "_TransmissionLinePad1";
  ValidationPad * pad1 =
    new ValidationPad(0,G4ThreeVector(dp_transmissionLinePad1Offset,0,0),
                      pad1Name,log_baseAlLayer,false,0,LM,
                      logicalLatticeContainer,borderContainer,pSurfChk);
  
  //Loop through the fundamental sub-volumes and push them back into the
  //fundamental subvolume list for the transmission line. We have to do this
  //here because the pads are composite volumes
  AddComplexGeometryPadSubVolumesToThisList(pad1);
  
  //Pad 2: rotate around Z axis by 180 degrees
  G4String pad2Name = pName + "_TransmissionLinePad2";
  G4RotationMatrix * pad2Rot = new G4RotationMatrix();
  pad2Rot->rotateZ(180*deg);
  ValidationPad * pad2 =
    new ValidationPad(pad2Rot,
                      G4ThreeVector(dp_transmissionLinePad2Offset,0,0),
                      pad2Name,log_baseAlLayer,false,0,LM,
                      logicalLatticeContainer,borderContainer,pSurfChk);
  
  //Loop through the fundamental sub-volumes and push them back into the
  //fundamental subvolume list for the transmission line. We have to do this
  //here because the pads are composite volumes
  AddComplexGeometryPadSubVolumesToThisList(pad2);
  
  //Now make the transmission line itself, in two parts: empty and conductor.
  //The logical mother volume of this is the base Ni layer.
  G4String tlNameEmpty = pName + "_TransmissionLineEmpty";
  G4String tlNameEmptySolid = tlNameEmpty + "_solid";
  G4String tlNameEmptyLog = tlNameEmpty + "_log";
  G4Box * solid_transmissionLineEmpty =
    new G4Box(tlNameEmptySolid,
              0.5 * (dp_transmissionLinePad2Offset -
                     dp_transmissionLinePad1Offset -
                     2 * dp_padEmptyPart2TrdZ - 2 * 0.5 * dp_padEmptyPart1DimX),
              0.5 * dp_transmissionLineCavityFullWidth,
              0.5 * dp_transmissionLineBaseLayerDimZ);
  
  G4LogicalVolume * log_transmissionLineEmpty
    = new G4LogicalVolume(solid_transmissionLineEmpty,air_mat,tlNameEmptyLog);  
  log_transmissionLineEmpty->SetVisAttributes(air_vis);
  
  G4VPhysicalVolume * phys_transmissionLineEmpty
    = new G4PVPlacement(0,G4ThreeVector(0,0,0),log_transmissionLineEmpty,
                        tlNameEmpty,log_baseAlLayer,false,0,true);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Vacuum",tlNameEmpty,phys_transmissionLineEmpty));

  
  G4String tlNameConductor = pName + "_TransmissionLineConductor";
  G4String tlNameConductorSolid = tlNameConductor + "_solid";
  G4String tlNameConductorLog = tlNameConductor + "_log";
  G4Box * solid_transmissionLineConductor =
    new G4Box(tlNameConductorSolid,
              0.5 * (dp_transmissionLinePad2Offset -
                     dp_transmissionLinePad1Offset -
                     2 * dp_padEmptyPart2TrdZ - 2 * 0.5 * dp_padEmptyPart1DimX),
              0.5 * dp_transmissionLineConductorWidth,
              0.5 * dp_transmissionLineBaseLayerDimZ);
  
  G4LogicalVolume * log_transmissionLineConductor
    = new G4LogicalVolume(solid_transmissionLineConductor,aluminum_mat,
                          tlNameConductorLog);
  log_transmissionLineConductor->SetVisAttributes(aluminum_vis);
  
  G4VPhysicalVolume * phys_transmissionLineConductor
    = new G4PVPlacement(0,G4ThreeVector(0,0,0),log_transmissionLineConductor,
                        tlNameConductor,log_transmissionLineEmpty,false,0,
                        true);
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Aluminum",tlNameConductor,phys_transmissionLineConductor));
  
  //Create dedicated lattice info for this TL piece
  G4LatticePhysical* AlPhysical_transmissionLineConductor
    = new G4LatticePhysical(AlLogical,dp_polycryElScatMFP_Al,
                            dp_scDelta0_Al, dp_scTeff_Al,
                            dp_scDn_Al, dp_scTauQPTrap_Al);
  AlPhysical_transmissionLineConductor->SetMillerOrientation(1,0,0);
  LM->RegisterLattice(phys_transmissionLineConductor,AlPhysical_transmissionLineConductor);  

  //Handle the intra-TL boundary setting:
  //1. pad 1 conductor to TL conductor (bidirectional)
  //2. pad 1 empty to TL empty (bidirectional)
  //3. pad 2 conductor to TL conductor (bidirectional)
  //4. pad 2 empty to TL empty (bidirectional)
  //5. TL conductor to TL empty
  //6. TL empty to TL conductor
  if (borderContainer.count("AlAl") == 0) {
    std::cout << "Uh oh. Trying to access borderContainer[AlAl] but it's not "
              << "there..." << std::endl;
  }
  G4CMPSurfaceProperty* AlAlBoundary = borderContainer["AlAl"];
  if (borderContainer.count("VacVac") == 0) {
    std::cout << "Uh oh. Trying to access borderContainer[VacVac] but it's not "
              << "there..." << std::endl;
  }
  G4CMPSurfaceProperty* VacVacBoundary = borderContainer["VacVac"];  
  if (borderContainer.count("AlVac") == 0) {
    std::cout << "Uh oh. Trying to access borderContainer[AlVac] but it's not "
              << "there..." << std::endl;
  }
  G4CMPSurfaceProperty* AlVacBoundary = borderContainer["AlVac"];
  G4String pad1TL_conductor_name = pName + "_Pad1TL_AlAl";
  G4String TLpad1_conductor_name = pName + "_TLPad1_AlAl";
  G4String pad1TL_empty_name = pName + "_Pad1TL_VacVac";
  G4String TLpad1_empty_name = pName + "_TLPad1_VacVac";
  G4String pad2TL_conductor_name = pName + "_Pad2TL_AlAl";
  G4String TLpad2_conductor_name = pName + "_TLPad2_AlAl";
  G4String pad2TL_empty_name = pName + "_Pad2TL_VacVac";
  G4String TLpad2_empty_name = pName + "_TLPad2_VacVac";
  G4String TLemptyTLconductor_name1 = pName + "_TLemptyTLconductor_AlVac";
  G4String TLemptyTLconductor_name2 = pName + "_TLemptyTLconductor_VacAl";
  
  //Find the volumes of interest
  std::map<std::string,G4VPhysicalVolume*> tempContainer;
  for (unsigned int iV = 0; iV < fFundamentalVolumeList.size(); ++iV) {
    if (std::get<1>(fFundamentalVolumeList[iV]).
        contains("TransmissionLinePad1_PadConductor")) {
      tempContainer.emplace("Pad1Conductor",
                            std::get<2>(fFundamentalVolumeList[iV]));
    }
    if (std::get<1>(fFundamentalVolumeList[iV]).
        contains("TransmissionLinePad1_PadEmpty")) {
      tempContainer.emplace("Pad1Empty",
                            std::get<2>(fFundamentalVolumeList[iV]));
    }    
    if (std::get<1>(fFundamentalVolumeList[iV]).
        contains("TransmissionLinePad2_PadConductor")) {
      tempContainer.emplace("Pad2Conductor",
                            std::get<2>(fFundamentalVolumeList[iV]));
    }
    if (std::get<1>(fFundamentalVolumeList[iV]).
        contains("TransmissionLinePad2_PadEmpty")) {
      tempContainer.emplace("Pad2Empty",
                            std::get<2>(fFundamentalVolumeList[iV]));
    }
  }
    

  
  new G4CMPLogicalBorderSurface(pad1TL_conductor_name,
                                phys_transmissionLineConductor,
                                tempContainer["Pad1Conductor"],AlAlBoundary);
  new G4CMPLogicalBorderSurface(TLpad1_conductor_name,
                                tempContainer["Pad1Conductor"],
                                phys_transmissionLineConductor,AlAlBoundary);
  
  new G4CMPLogicalBorderSurface(pad1TL_empty_name,
                                phys_transmissionLineEmpty,
                                tempContainer["Pad1Empty"],VacVacBoundary);
  new G4CMPLogicalBorderSurface(TLpad1_empty_name,
                                tempContainer["Pad1Empty"],
                                phys_transmissionLineEmpty,VacVacBoundary);

  new G4CMPLogicalBorderSurface(pad2TL_conductor_name,
                                phys_transmissionLineConductor,
                                tempContainer["Pad2Conductor"],AlAlBoundary);
  new G4CMPLogicalBorderSurface(TLpad2_conductor_name,
                                tempContainer["Pad2Conductor"],
                                phys_transmissionLineConductor,AlAlBoundary);
  
  new G4CMPLogicalBorderSurface(pad2TL_empty_name,
                                phys_transmissionLineEmpty,
                                tempContainer["Pad2Empty"],VacVacBoundary);
  new G4CMPLogicalBorderSurface(TLpad2_empty_name,
                                tempContainer["Pad2Empty"],
                                phys_transmissionLineEmpty,VacVacBoundary);

  new G4CMPLogicalBorderSurface(TLemptyTLconductor_name1,
                                phys_transmissionLineConductor,
                                phys_transmissionLineEmpty,AlVacBoundary);
  new G4CMPLogicalBorderSurface(TLemptyTLconductor_name2,
                                phys_transmissionLineEmpty,
                                phys_transmissionLineConductor,AlVacBoundary);

  


  ///////////////////////////////////////////
  // Output logical/physical volume selection
  //-----------------------------------------

  fLog_output = log_baseAlLayer;
  fPhys_output = phys_baseAlLayer;


}

  

//Used for building the base layer for the Aluminum. Can't be a pure rectangle
//because that rectangle will overlap with the resonators, which we can't
//have... Ripping some of this from the code in the pad file and the
//transmission line code, so be careful if you end up having to change the
//dimensions of things...
G4UnionSolid * ValidationTransmissionLine::
CreatePieceBasedAlLayer(G4String nameSolid) {

  //For the pad, we start with a volume of air with the same dimensions as the
  //aluminum substrate. Then we add the aluminum pad within.
  G4Box * solid_padEmptyPart1
    = new G4Box("baseAlLayerEmptyPadPart1Solid",0.5 * dp_padEmptyPart1DimX,
                0.5 * dp_padEmptyPart1DimY,0.5 * dp_padEmptyPart1DimZ);
  G4Trd * solid_padEmptyPart2
    = new G4Trd("baseAlLayerEmptyPadPart2Solid", 0.5 * dp_padEmptyPart2TrdX1,
                0.5 * dp_padEmptyPart2TrdX2,0.5 * dp_padEmptyPart2TrdY1,
                0.5 * dp_padEmptyPart2TrdY2,0.5 * dp_padEmptyPart2TrdZ);
  
  //We need to rotate the part 2 so that it can be aligned and placed next to
  //part 1
  G4RotationMatrix * rotEmptyPart2 = new G4RotationMatrix();
  rotEmptyPart2->rotateX(90.*deg);
  rotEmptyPart2->rotateY(-90.*deg);

  G4ThreeVector transPart2EmptyWrtPart1Empty(dp_padEmptyPart1DimX/2.0
                                             + dp_padEmptyPart2TrdZ/2.0,0,0);
  G4UnionSolid * solid_padEmpty =
    new G4UnionSolid("BaseAlLayerEmptyPadSolid",solid_padEmptyPart1,
                     solid_padEmptyPart2,rotEmptyPart2,
                     transPart2EmptyWrtPart1Empty);
  
  G4Box * solid_transmissionLineEmpty
    = new G4Box("BaseAlLayerTransmissionLineSolid",
                0.5 * (dp_transmissionLinePad2Offset -
                       dp_transmissionLinePad1Offset -
                       2 * dp_padEmptyPart2TrdZ -
                       2 * 0.5 * dp_padEmptyPart1DimX),
                0.5 * dp_transmissionLineCavityFullWidth,
                0.5 * dp_transmissionLineBaseLayerDimZ);
  
  G4UnionSolid * solid_merger1 =
    new G4UnionSolid("Merger1",solid_transmissionLineEmpty,solid_padEmpty,0,
                     G4ThreeVector(dp_transmissionLinePad1Offset,0,0));

  G4RotationMatrix * pad2Rot = new G4RotationMatrix();
  pad2Rot->rotateZ(180*deg);
  G4UnionSolid * solid_baseAlLayer =
    new G4UnionSolid(nameSolid,solid_merger1,solid_padEmpty,pad2Rot,
                     G4ThreeVector(dp_transmissionLinePad2Offset,0,0));

  return solid_baseAlLayer;
  
}

std::vector<std::tuple<std::string,G4String,G4VPhysicalVolume*> >
ValidationTransmissionLine::GetListOfAllFundamentalSubVolumes() {
  return fFundamentalVolumeList;
}

void ValidationTransmissionLine::
AddComplexGeometryPadSubVolumesToThisList(ValidationPad * pad) {
  for (unsigned int iSubVol = 0; iSubVol < pad->GetListOfAllFundamentalSubVolumes().size(); ++iSubVol) {
    std::tuple<std::string,G4String,G4VPhysicalVolume*>
      theTuple(std::get<0>(pad->GetListOfAllFundamentalSubVolumes()[iSubVol]),
               std::get<1>(pad->GetListOfAllFundamentalSubVolumes()[iSubVol]),
               std::get<2>(pad->GetListOfAllFundamentalSubVolumes()[iSubVol]));
    fFundamentalVolumeList.push_back(theTuple);
  }
}

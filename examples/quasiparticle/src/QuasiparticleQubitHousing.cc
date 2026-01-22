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
// 20260109  M. Kelsey -- G4CMP-569: Removed unused local variables

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
#include "QuasiparticleQubitHousing.hh"
#include "QuasiparticleDetectorParameters.hh"

using namespace QuasiparticleDetectorParameters;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Primary Constructor
QuasiparticleQubitHousing::
QuasiparticleQubitHousing(G4RotationMatrix* pRot, const G4ThreeVector& tLate,
                          const G4String & pName,
                          G4LogicalVolume* pMotherLogical, G4bool pMany,
                          G4int pCopyNo,
                          G4bool pSurfChk) {
  //Here, use the inputs to this to set up the geometry and fill out the
  //PVPlacement data member, which is the real output from this class (and
  //which we'll access in our main detector construction file.)
  ConstructQubitHousing(pRot,tLate,
                        pName,
                        pMotherLogical,
                        pMany,
                        pCopyNo,
                        pSurfChk);
  
  
}

// Default Constructor
QuasiparticleQubitHousing::QuasiparticleQubitHousing() {
}

// Destructor
QuasiparticleQubitHousing::~QuasiparticleQubitHousing() {
}

//Moving implementation down here so it's not in the constructor
void QuasiparticleQubitHousing::
ConstructQubitHousing(G4RotationMatrix * pRot,const G4ThreeVector & tLate,
                      const G4String & pName,G4LogicalVolume * pMotherLogical,
                      G4bool pMany,G4int pCopyNo,G4bool pSurfChk) {

  //Start with some preliminaries - NIST manager
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* copper_mat = nist->FindOrBuildMaterial("G4_Cu");
  
  //First, start with the geometrical objects / solid volumes of relevance

  //Start with solid block of copper, then carve out the stuff where there is
  //empty space in the housing. The radial cutouts we are treating
  //as rectangular prisms for now (close enough to reality)
  G4Box * solid_qubitHousingBlock = new G4Box("QubitHousingBlock_solid",
                                              0.5*dp_housingDimX,
                                              0.5*dp_housingDimY,
                                              0.5*dp_housingDimZ);
  
  G4Box * solid_centralCutout = new G4Box("CentralCutout_solid",
                                          0.5*dp_housingCentralCutoutDimX,
                                          0.5*dp_housingCentralCutoutDimY,
                                          0.5*dp_housingCentralCutoutDimZ);
  
  G4Box * solid_radialCutout = new G4Box("RadialCutout_solid",
                                         0.5 * dp_housingRadialCutoutDimX,
                                         0.5 * dp_housingRadialCutoutDimY,
                                         0.5 * dp_housingRadialCutoutDimZ);
  
  

  //-----------------------------------------------------------------------
  //Now take the radial cutouts and move them to the right position relative
  //to the primary cutout

  //First of eight
  G4RotationMatrix * rot1 = new G4RotationMatrix();
  rot1->rotateZ(45.*deg);
  G4ThreeVector trans1(0.5*dp_housingCentralCutoutDimX, //Top right corner
                       0.5*dp_housingCentralCutoutDimY, //Top right corner
                       0.5*(dp_housingCentralCutoutDimZ
                            - dp_housingRadialCutoutDimZ));//Top of chip  
  G4UnionSolid * solid_partialCutout1 =
    new G4UnionSolid("PartialCutout1_solid", solid_centralCutout,
                     solid_radialCutout, rot1, trans1);
  
  //Second of eight
  G4RotationMatrix * rot2 = new G4RotationMatrix();
  rot2->rotateZ(-45.*deg);
  G4ThreeVector trans2(0.5*dp_housingCentralCutoutDimX, //Bottom right corner
		       -0.5*dp_housingCentralCutoutDimY, //Bottom right corner
		       0.5*(dp_housingCentralCutoutDimZ
                    - dp_housingRadialCutoutDimZ));//Top of chip  
  G4UnionSolid * solid_partialCutout2 =
    new G4UnionSolid("PartialCutout2_solid",solid_partialCutout1,
                     solid_radialCutout,rot2,trans2);

  //Third of eight
  G4RotationMatrix * rot3 = new G4RotationMatrix();
  rot3->rotateZ(-45.*deg);
  G4ThreeVector trans3(-0.5*dp_housingCentralCutoutDimX, //Top left corner
                       0.5*dp_housingCentralCutoutDimY, //Top left corner
                       0.5*(dp_housingCentralCutoutDimZ
                            - dp_housingRadialCutoutDimZ));//Move to top of chip  
  G4UnionSolid * solid_partialCutout3 =
    new G4UnionSolid("PartialCutout3_solid",solid_partialCutout2,
                     solid_radialCutout,rot3,trans3);
  
  //Fourth of eight
  G4RotationMatrix * rot4 = new G4RotationMatrix();
  rot4->rotateZ(45.*deg);
  G4ThreeVector trans4(-0.5*dp_housingCentralCutoutDimX, //Bottom left corner
                       -0.5*dp_housingCentralCutoutDimY, //Bottom left corner
                       0.5*(dp_housingCentralCutoutDimZ
                            - dp_housingRadialCutoutDimZ));//Top of chip  
  G4UnionSolid * solid_partialCutout4 =
    new G4UnionSolid("PartialCutout4_solid",solid_partialCutout3,
                     solid_radialCutout,rot4,trans4);


  //Fifth of eight
  G4RotationMatrix * rot5 = new G4RotationMatrix();
  rot5->rotateZ(90.*deg);
  G4ThreeVector trans5(0.5*dp_housingCentralCutoutDimX, //Move to top in X
                       0,0.5*(dp_housingCentralCutoutDimZ -
                              dp_housingRadialCutoutDimZ));//Top of chip  
  G4UnionSolid * solid_partialCutout5 =
    new G4UnionSolid("PartialCutout5_solid",solid_partialCutout4,
                     solid_radialCutout,rot5,trans5);
  
  //Sixth of eight
  G4RotationMatrix * rot6 = new G4RotationMatrix();
  rot6->rotateZ(90.*deg);
  G4ThreeVector trans6(-0.5*dp_housingCentralCutoutDimX, //Bottom in X
                       0,0.5*(dp_housingCentralCutoutDimZ -
                              dp_housingRadialCutoutDimZ));//Top of chip  
  G4UnionSolid * solid_partialCutout6 =
    new G4UnionSolid("PartialCutout6_solid",solid_partialCutout5,
                     solid_radialCutout,rot6,trans6);
 
  //Seventh of eight
  G4RotationMatrix * rot7 = new G4RotationMatrix();
  rot7->rotateZ(0.*deg);
  G4ThreeVector trans7(0,0.5*dp_housingCentralCutoutDimY, //Move to bottom (Y)
                       0.5*(dp_housingCentralCutoutDimZ -
                            dp_housingRadialCutoutDimZ)); //Top of chip  
  G4UnionSolid * solid_partialCutout7 =
    new G4UnionSolid("PartialCutout7_solid",solid_partialCutout6,
                     solid_radialCutout,rot7,trans7);
  
  //Eigth of eight
  G4RotationMatrix * rot8 = new G4RotationMatrix();
  rot8->rotateZ(0.*deg);
  G4ThreeVector trans8(0,-0.5*dp_housingCentralCutoutDimY, //Move to top (Y)
                       0.5*(dp_housingCentralCutoutDimZ -
                            dp_housingRadialCutoutDimZ));//Move to top of chip  
  G4UnionSolid * solid_totalCutout =
    new G4UnionSolid("TotalCutout_solid",solid_partialCutout7,
                     solid_radialCutout,rot8,trans8);


  //-----------------------------------------------------------------------
  //Now take the total cutout and subtract it from the main block
  //Modifying by one micron so that we're unambiguously removing the top
  //surface of housing
  G4ThreeVector
    cutoutTranslation(0.0,0.0,(dp_housingDimZ-
                               dp_housingCentralCutoutDimZ)/2.0+dp_eps);
  
  G4SubtractionSolid * solid_qubitHousingNoAddIn =
    new G4SubtractionSolid("QubitHousingNoAddIn_solid",solid_qubitHousingBlock,
                           solid_totalCutout,0,cutoutTranslation);
  

  //...
  //Do any boolean combination of multiple geometric objects here, I think.
  //...

  
  //-----------------------------------------------------------------------
  //Need to add back in corners under the cutouts. Need to do some math here to get this right...
  G4Box * solid_radialCornerAddIn =
    new G4Box("RadialCornerAddIn_solid",
              0.5 * dp_housingRadialCutoutDimX,0.5 * dp_housingRadialCutoutDimX,
              (0.5 * (dp_housingCentralCutoutDimZ -
                      dp_housingRadialCutoutDimZ)));

  //First corner
  G4RotationMatrix * rot9 = new G4RotationMatrix();
  rot9->rotateZ(45.*deg);  
  G4ThreeVector addInTranslation1(0.5*dp_housingCentralCutoutDimX,
                                  0.5*dp_housingCentralCutoutDimY,
                                  -(dp_housingRadialCutoutDimZ)/2.0);
  G4UnionSolid * solid_qubitHousingAddIn1 =
    new G4UnionSolid("QubitHousingAddIn1_solid",solid_qubitHousingNoAddIn,
                     solid_radialCornerAddIn,rot9,addInTranslation1 +
                     cutoutTranslation);

  //Second corner
  G4ThreeVector addInTranslation2(0.5*dp_housingCentralCutoutDimX,
                                  -0.5*dp_housingCentralCutoutDimY,
                                  -(dp_housingRadialCutoutDimZ/2.0));
  G4UnionSolid * solid_qubitHousingAddIn2 =
    new G4UnionSolid("QubitHousingAddIn2_solid",solid_qubitHousingAddIn1,
                     solid_radialCornerAddIn, rot9, addInTranslation2 +
                     cutoutTranslation);


  //Third corner
  G4ThreeVector addInTranslation3(-0.5*dp_housingCentralCutoutDimX,
                                  -0.5*dp_housingCentralCutoutDimY,
                                  -(dp_housingRadialCutoutDimZ/2.0));
  G4UnionSolid * solid_qubitHousingAddIn3 =
    new G4UnionSolid("QubitHousingAddIn3_solid",solid_qubitHousingAddIn2,
                     solid_radialCornerAddIn,rot9,addInTranslation3 +
                     cutoutTranslation);
  
  //Fourth corner
  G4ThreeVector addInTranslation4(-0.5*dp_housingCentralCutoutDimX,
                                  0.5*dp_housingCentralCutoutDimY,
                                  -(dp_housingRadialCutoutDimZ/2.0));
  G4UnionSolid * solid_qubitHousing =
    new G4UnionSolid("QubitHousing_solid",solid_qubitHousingAddIn3,
                     solid_radialCornerAddIn,rot9,addInTranslation4 +
                     cutoutTranslation);

 
  ///////////////////////////////////////////
  // Logical and Physical Volume Creation
  //-----------------------------------------
  // Done after one solid object subsuming all of the geometry has been made.
  
  //Now attribute a physical material to the housing
  G4LogicalVolume * log_QubitHousing =
    new G4LogicalVolume(solid_qubitHousing, copper_mat, pName);
  
  //Now, create a physical volume and G4PVPlacement for storing as the final
  //output
  G4VPhysicalVolume* phys_QubitHousing =
    new G4PVPlacement(pRot, tLate, log_QubitHousing, pName, pMotherLogical,
                      pMany, pCopyNo, pSurfChk);




  G4VisAttributes* simpleBoxVisAtt =
    new G4VisAttributes(G4Colour(1.0,0.647,0.0,0.9));
  simpleBoxVisAtt->SetVisibility(true);
  log_QubitHousing->SetVisAttributes(simpleBoxVisAtt);
  
  // Lastly, make the logical volume and physical volume accessible data members
  fLog_output = log_QubitHousing;
  fPhys_output = phys_QubitHousing; 
}

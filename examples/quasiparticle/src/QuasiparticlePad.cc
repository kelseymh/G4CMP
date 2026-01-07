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
#include "G4LatticePhysical.hh"
#include "G4CMPLogicalBorderSurface.hh"

//Includes (specific to this project)
#include "QuasiparticlePad.hh"
#include "QuasiparticleDetectorParameters.hh"

using namespace QuasiparticleDetectorParameters;

// Primary Constructor
QuasiparticlePad::
QuasiparticlePad(G4RotationMatrix * pRot,const G4ThreeVector & tLate,
		 const G4String & pName,G4LogicalVolume * pMotherLogical,
		 G4bool pMany,G4int pCopyNo,G4LatticeManager * LM,
		 std::map<std::string,G4LatticeLogical*> logicalLatticeContainer,
		 std::map<std::string,G4CMPSurfaceProperty*> borderContainer,
		 G4bool pSurfChk) {
  //Here, use the inputs to this to set up the geometry and fill out the
  //PVPlacement data member, which is the real output from this class (and
  //which we'll access in our main detector construction
  //file.)
  ConstructPad(pRot,tLate,pName,pMotherLogical,pMany,pCopyNo,LM,
	       logicalLatticeContainer,borderContainer,pSurfChk);
}

// Default Constructor
QuasiparticlePad::QuasiparticlePad() {  
}

// Destructor
QuasiparticlePad::~QuasiparticlePad() {
}

//Moving implementation down here so it's not in the constructor
void QuasiparticlePad::
ConstructPad(G4RotationMatrix * pRot,const G4ThreeVector & tLate,
	     const G4String & pName,G4LogicalVolume * pMotherLogical,
	     G4bool pMany,G4int pCopyNo,G4LatticeManager * LM,
	     std::map<std::string,G4LatticeLogical*> logicalLatticeContainer,
	     std::map<std::string,G4CMPSurfaceProperty*> borderContainer,
	     G4bool pSurfChk) {

  //Start with some preliminaries - NIST manager
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* aluminum_mat = nist->FindOrBuildMaterial("G4_Al");
  G4Material* air_mat = nist->FindOrBuildMaterial("G4_AIR");

  //Set up lattice information
  G4LatticeLogical* AlLogical = logicalLatticeContainer["Aluminum"];
  
  //Set up the aluminum visualization
  G4VisAttributes* aluminum_vis= new G4VisAttributes(G4Colour(0.0,1.0,1.0,0.5));
  aluminum_vis->SetVisibility(true);
  G4VisAttributes* air_vis= new G4VisAttributes(G4Colour(0.5,0.5,0.5,0.5));
  air_vis->SetVisibility(true);
  

  //--------------------------------------------------------------------------
  //For the pad, we start with a volume of air with the same dimensions as the
  //aluminum substrate. Then we add the aluminum pad within.
  G4String padEmptyPart1Name = pName + "_PadEmptyPart1";
  G4String padEmptyPart1NameLog = pName + "_PadEmptyPart1_log";
  G4String padEmptyPart1NameSolid = pName + "_PadEmptyPart1_solid";
  G4String padEmptyPart2Name = pName + "_PadEmptyPart2";
  G4String padEmptyPart2NameLog = pName + "_PadEmptyPart2_log";
  G4String padEmptyPart2NameSolid = pName + "_PadEmptyPart2_solid";
  G4Box * solid_padEmptyPart1 =
    new G4Box(padEmptyPart1NameSolid,0.5 * dp_padEmptyPart1DimX,
	      0.5 * dp_padEmptyPart1DimY,0.5 * dp_padEmptyPart1DimZ);
  G4Trd * solid_padEmptyPart2 =
    new G4Trd(padEmptyPart2NameSolid,0.5 * dp_padEmptyPart2TrdX1,
	      0.5 * dp_padEmptyPart2TrdX2,0.5 * dp_padEmptyPart2TrdY1,
	      0.5 * dp_padEmptyPart2TrdY2,0.5 * dp_padEmptyPart2TrdZ);

  //We need to rotate the part 2 so that it can be aligned and placed next to
  //part 1
  G4RotationMatrix * rotEmptyPart2 = new G4RotationMatrix();
  rotEmptyPart2->rotateX(90.*deg);
  rotEmptyPart2->rotateY(-90.*deg);

  G4ThreeVector transPart2EmptyWrtPart1Empty(dp_padEmptyPart1DimX/2.0
					     + dp_padEmptyPart2TrdZ/2.0,0,0);

  G4String padEmptyName = pName + "_PadEmpty";
  G4String padEmptyNameLog = pName + "_PadEmpty_log";
  G4String padEmptyNameSolid = pName + "_PadEmpty_solid";
  G4UnionSolid * solid_padEmpty =
    new G4UnionSolid(padEmptyNameSolid,solid_padEmptyPart1,solid_padEmptyPart2,
		     rotEmptyPart2,transPart2EmptyWrtPart1Empty);

  //Now we make a smaller pad that is internal to the larger one. This is
  //actually aluminum and is where the electrical contact is actually made
  G4String padConductorPart1Name = pName + "_PadConductorPart1";
  G4String padConductorPart1NameLog = pName + "_PadConductorPart1_log";
  G4String padConductorPart1NameSolid = pName + "_PadConductorPart1_solid";
  G4String padConductorPart2Name = pName + "_PadConductorPart2";
  G4String padConductorPart2NameLog = pName + "_PadConductorPart2_log";
  G4String padConductorPart2NameSolid = pName + "_PadConductorPart2_solid";
  G4Box * solid_padConductorPart1 =
    new G4Box(padConductorPart1NameSolid,0.5 * dp_padPart1DimX,
	      0.5 * dp_padPart1DimY,0.5 * dp_padPart1DimZ);
  G4Trd * solid_padConductorPart2 =
    new G4Trd(padConductorPart2NameSolid,0.5 * dp_padPart2TrdX1,
	      0.5 * dp_padPart2TrdX2,0.5 * dp_padPart2TrdY1,
	      0.5 * dp_padPart2TrdY2,0.5 * dp_padPart2TrdZ);

  
  //We need to rotate the part 2 so that it can be aligned and placed next to
  //part 1
  G4RotationMatrix * rotPart2 = new G4RotationMatrix();
  rotPart2->rotateX(90.*deg);
  rotPart2->rotateY(-90.*deg);

  G4ThreeVector transPart2WrtPart1(dp_padPart1DimX/2.0
				   + dp_padPart2TrdZ/2.0,0,0);

  G4String padConductorName = pName + "_PadConductor";
  G4String padConductorNameLog = pName + "_PadConductor_log";
  G4String padConductorNameSolid = pName + "_PadConductor_solid";
  G4UnionSolid * solid_padConductor =
    new G4UnionSolid(padConductorNameSolid,solid_padConductorPart1,
		     solid_padConductorPart2,rotPart2,transPart2WrtPart1);
  
  //We now create a shift the non-empty pad and place it inside:
  G4ThreeVector transPadWrtEmptyPad(dp_padPart2InternalShiftX,0,0);
 
  ///////////////////////////////////////////
  // Logical and Physical Volume Creation
  //-----------------------------------------
  // Done after one solid object subsuming all of the geometry has been made.

  //Now attribute a physical material to the housing
  G4LogicalVolume * log_padEmpty =
    new G4LogicalVolume(solid_padEmpty,air_mat,padEmptyNameLog);
  log_padEmpty->SetVisAttributes(air_vis);
  
  G4LogicalVolume * log_padConductor =
    new G4LogicalVolume(solid_padConductor,aluminum_mat,padConductorNameLog);
  log_padConductor->SetVisAttributes(aluminum_vis);
  
  //Now, create a physical volume and G4PVPlacement for storing as the final
  //output. Also, create a dedicated lattice for this.
  G4VPhysicalVolume* phys_padConductor =
    new G4PVPlacement(0,transPadWrtEmptyPad,log_padConductor,padConductorName,
		      log_padEmpty,false,0,true);
  G4LatticePhysical* AlPhysical_padConductor =
    new G4LatticePhysical(AlLogical,dp_polycryElScatMFP_Al,
			  dp_scDelta0_Al,dp_scTeff_Al,dp_scDn_Al,
			  dp_scTauQPTrap_Al);
  AlPhysical_padConductor->SetMillerOrientation(1,0,0);
  LM->RegisterLattice(phys_padConductor,AlPhysical_padConductor);  
  G4VPhysicalVolume* phys_padEmpty =
    new G4PVPlacement(pRot,tLate,log_padEmpty,padEmptyName,pMotherLogical,
		      pMany,pCopyNo,pSurfChk);


  //Establish intra-pad boundaries
  if (borderContainer.count("AlVac") == 0){
    std::cout << "Uh oh. Trying to access borderContainer[AlVac] but it's not "
	      << "there..." << std::endl;
  }
  G4CMPSurfaceProperty* AlVacBoundary = borderContainer["AlVac"];
  G4String boundaryName1 = pName + "_AlVac";
  G4String boundaryName2 = pName + "_VacAl";
  new G4CMPLogicalBorderSurface(boundaryName1, phys_padConductor, phys_padEmpty,
				AlVacBoundary);
  new G4CMPLogicalBorderSurface(boundaryName2, phys_padEmpty, phys_padConductor,
				AlVacBoundary);

  
  //Push these back into the fundamental volume list. These will also serve as
  //the things we need to establish boundaries for
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Vacuum",padEmptyName,phys_padEmpty));
  fFundamentalVolumeList.push_back(std::tuple<std::string,G4String,G4VPhysicalVolume*>("Aluminum",padConductorName,phys_padConductor));
  
  // Lastly, make the logical volume and physical volume of the top/mother object accessible data members
  fLog_output = log_padEmpty;
  fPhys_output = phys_padEmpty;
}

std::vector<std::tuple<std::string,G4String,G4VPhysicalVolume*> >
QuasiparticlePad::GetListOfAllFundamentalSubVolumes() {
  return fFundamentalVolumeList;
}

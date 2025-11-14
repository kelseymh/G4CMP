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
/// \file ValidationPad.hh
/// \brief Definition of the class defining how a transmission line pad is constructed on a qubit chip
///        

#ifndef ValidationPad_h
#define ValidationPad_h 1


//#include "G4PVPlacement.hh"
#include "ValidationDetectorParameters.hh"
#include "G4LatticeManager.hh"
#include "G4LatticeLogical.hh"
#include "G4CMPSurfaceProperty.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4PVPlacement;

/// Detector construction class to define materials and geometry.
class ValidationPad
{
  public:
    ValidationPad();
    ~ValidationPad();

    //This is the constructor that should be used in general. It does not have the current logical
    //volume included because that will be defined IN the Qubit housing implementation. All we need is
    //a set of info that is external to this, which should be self-contained.
    ValidationPad(G4RotationMatrix * pRot,
		     const G4ThreeVector & tLate,
		     const G4String & pName,
		     G4LogicalVolume * pMotherLogical,
		     G4bool pMany,
		     G4int pCopyNo,
		     G4LatticeManager * LM,
		     std::map<std::string,G4LatticeLogical*> logicalLatticeContainer,
		     std::map<std::string,G4CMPSurfaceProperty*> borderContainer,
		     G4bool pSurfChk=false);

  
    //Access functions
    G4VPhysicalVolume * GetPhysicalVolume() { return fPhys_output; }
    G4LogicalVolume * GetLogicalVolume() { return fLog_output; }
  
    //Misc
    void ConstructPad(G4RotationMatrix * pRot,
		      const G4ThreeVector & tLate,
		      const G4String & pName,
		      G4LogicalVolume * pMotherLogical,
		      G4bool pMany,
		      G4int pCopyNo,
		      G4LatticeManager * LM,
		      std::map<std::string,G4LatticeLogical*> logicalLatticeContainer,
		      std::map<std::string,G4CMPSurfaceProperty*> borderContainer,
		      G4bool pSurfChk=false);
  
    std::vector<std::tuple<std::string,G4String,G4VPhysicalVolume*> > GetListOfAllFundamentalSubVolumes();
  
  protected:

  private:

    //The final G4PVPlacement
  G4LogicalVolume * fLog_output;
  G4VPhysicalVolume * fPhys_output;
  std::vector<std::tuple<std::string,G4String,G4VPhysicalVolume*> > fFundamentalVolumeList; //List of all fundamental sub-volumes in the pad
  
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

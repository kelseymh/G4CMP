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
/// \file QuasiparticleQubitHousing.hh
/// \brief Definition of the class defining shield type A for the
///        QUIET fridge

#ifndef QuasiparticleQubitHousing_h
#define QuasiparticleQubitHousing_h 1


//#include "G4PVPlacement.hh"
#include "QuasiparticleDetectorParameters.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4PVPlacement;

/// Detector construction class to define materials and geometry.
class QuasiparticleQubitHousing
{
  public:
    QuasiparticleQubitHousing();
    ~QuasiparticleQubitHousing();

    //This is the constructor that should be used in general. It does not have
    //the current logical volume included because that will be defined IN the
    //Qubit housing implementation. All we need is a set of info that is
    //external to this, which should be self-contained.
    QuasiparticleQubitHousing(G4RotationMatrix* pRot,
                              const G4ThreeVector& tLate,
                              const G4String& pName,
                              G4LogicalVolume* pMotherLogical,
                              G4bool pMany,
                              G4int pCopyNo,
                              G4bool pSurfChk=false);
  
    //Access functions
    G4VPhysicalVolume * GetPhysicalVolume(){ return fPhys_output; }
    G4LogicalVolume * GetLogicalVolume(){ return fLog_output; }
  
    //Misc
    void ConstructQubitHousing(G4RotationMatrix* pRot,
                               const G4ThreeVector& tLate,
                               const G4String& pName,
                               G4LogicalVolume* pMotherLogical,
                               G4bool pMany,
                               G4int pCopyNo,
                               G4bool pSurfChk=false);
  
  
  protected:

  private:

    //The final G4PVPlacement
    G4LogicalVolume* fLog_output;
    G4VPhysicalVolume* fPhys_output;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

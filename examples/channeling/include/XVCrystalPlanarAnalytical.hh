/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id: 2ef25be581a1c45efbae5deab6c3716825359f55 $
//
#ifndef XVCrystalPlanarAnalytical_h
#define XVCrystalPlanarAnalytical_h

#include "XVCrystalCharacteristic.hh"

class XVCrystalPlanarAnalytical:public XVCrystalCharacteristic {

private:
    G4int fNumberOfPlanes;

public:
    //set function
    void SetNumberOfPlanes(G4int);

    //retrieval function
    G4int GetNumberOfPlanes();

    virtual G4double ComputeECForSinglePlane(G4double,XPhysicalLattice*) = 0; // G4double = position in the channel

    //virtual function of XVCrystalCharacteristic
    G4ThreeVector ComputeEC(G4ThreeVector,XPhysicalLattice*);
    G4ThreeVector ComputeECFromVector(G4ThreeVector);
    G4ThreeVector ComputePositionInUnitCell(G4ThreeVector,XPhysicalLattice*);//G4double = position in the channel; G4double& = interplanar distance
    
    virtual G4double ComputeMaximum(XPhysicalLattice*);
    virtual G4double ComputeMinimum(XPhysicalLattice*);

    virtual void PrintOnFile(const G4String&,XPhysicalLattice*,G4double);
    
    void InitializeVector();
    //Contructors
    XVCrystalPlanarAnalytical();
    ~XVCrystalPlanarAnalytical();
};

#endif

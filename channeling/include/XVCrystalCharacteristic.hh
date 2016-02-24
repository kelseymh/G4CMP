/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id: 13851f75247fbf80f644b7ecfe1c2b6446a6b6ac $
//
#ifndef XVCrystalCharacteristic_h
#define XVCrystalCharacteristic_h

#include "XLatticeManager3.hh"
#include "G4Track.hh"
#include "G4PhysicsVector.hh"

class XVCrystalCharacteristic {

private:
    XLatticeManager3* fLatticeManager;

protected:
    G4double fMaximum;
    G4double fMinimum;
    XPhysicalLattice *fPhysicalLattice;
    G4PhysicsVector* fVectorEC;
    
public:
    //retrieval functions
    XPhysicalLattice* GetXPhysicalLattice(G4VPhysicalVolume*);
    XUnitCell* GetXUnitCell(G4VPhysicalVolume*);
    XLogicalLattice* GetLogicalLattice(G4VPhysicalVolume*);
    void InitializePhysicalLattice(XPhysicalLattice*);
    
    //virtual function to compute value starting from the point in the xtal reference frame and the physical volume of the xtal
    G4ThreeVector GetEC(G4ThreeVector,XPhysicalLattice*);
    virtual G4ThreeVector ComputeEC(G4ThreeVector,XPhysicalLattice*) = 0;
    virtual G4ThreeVector ComputeECFromVector(G4ThreeVector) = 0;
    virtual G4ThreeVector ComputePositionInUnitCell(G4ThreeVector,XPhysicalLattice*);
    
    virtual G4double ComputeTFScreeningRadius(XPhysicalLattice*);

    virtual G4double GetMaximum(XPhysicalLattice*);
    virtual G4double GetMinimum(XPhysicalLattice*);

    virtual G4double ComputeMaximum(XPhysicalLattice*);
    virtual G4double ComputeMinimum(XPhysicalLattice*);
    
    virtual void PrintOnFile(const G4String&,XPhysicalLattice*,G4double = 1) = 0;
    G4bool IsInitialized(XPhysicalLattice*);
    virtual void InitializeVector() = 0;
    
    //Contructors
    XVCrystalCharacteristic();
    ~XVCrystalCharacteristic();
};

#endif

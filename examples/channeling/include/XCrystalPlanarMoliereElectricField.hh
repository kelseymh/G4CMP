/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id: 3da6feea72e9277d5ce84781485ac7f98673a0a4 $
//
#ifndef XCrystalPlanarMoliereElectricField_h
#define XCrystalPlanarMoliereElectricField_h

#include "XVCrystalPlanarAnalytical.hh"

class XCrystalPlanarMoliereElectricField:public XVCrystalPlanarAnalytical {

private:
    G4double fAlfa[3];
    G4double fBeta[3];

public:
    G4double ComputeECForSinglePlane(G4double,XPhysicalLattice*);
    
    //Contructors
    XCrystalPlanarMoliereElectricField();
    ~XCrystalPlanarMoliereElectricField();
};

#endif

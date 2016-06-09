/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id: e2b75b494b983130d0f571a7abd3b5c0f4ff085e $
//
#ifndef XCrystalPlanarMolierePotential_h
#define XCrystalPlanarMolierePotential_h

#include "XVCrystalPlanarAnalytical.hh"

class XCrystalPlanarMolierePotential:public XVCrystalPlanarAnalytical {

private:
    G4double fAlfa[3];
    G4double fBeta[3];

public:
    G4double ComputeECForSinglePlane(G4double,XPhysicalLattice*);
    
    G4double ComputeMaximum(XPhysicalLattice*);
    G4double ComputeMinimum(XPhysicalLattice*);

    //Contructors
    XCrystalPlanarMolierePotential();
    ~XCrystalPlanarMolierePotential();
};

#endif

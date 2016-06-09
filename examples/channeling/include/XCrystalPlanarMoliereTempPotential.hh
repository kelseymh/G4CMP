/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id: 5c22785a0e0d30faaff99c11741501599d2c6aa3 $
//
#ifndef XCrystalPlanarMoliereTempPotential_h
#define XCrystalPlanarMoliereTempPotential_h

#include "XVCrystalPlanarAnalytical.hh"

class XCrystalPlanarMoliereTempPotential:public XVCrystalPlanarAnalytical {

private:
    G4double fAlfa[3];
    G4double fBeta[3];
    G4double fGamma[3];

public:
    G4double ComputeECForSinglePlane(G4double,XPhysicalLattice*);
    
    G4double ComputeMaximum(XPhysicalLattice*);
    G4double ComputeMinimum(XPhysicalLattice*);

    //Contructors
    XCrystalPlanarMoliereTempPotential();
    ~XCrystalPlanarMoliereTempPotential();
};

#endif

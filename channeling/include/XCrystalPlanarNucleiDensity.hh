/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id: 437fcbc996bace770da7fe579d1c1736d8a3bd14 $
//
#ifndef XCrystalPlanarNucleiDensity_h
#define XCrystalPlanarNucleiDensity_h

#include "XVCrystalPlanarAnalytical.hh"

class XCrystalPlanarNucleiDensity:public XVCrystalPlanarAnalytical {

private:

public:
    G4double ComputeECForSinglePlane(G4double,XPhysicalLattice*);
    
    //Contructors
    XCrystalPlanarNucleiDensity();
    ~XCrystalPlanarNucleiDensity();
};

#endif

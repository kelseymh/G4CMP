/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id: dd717f046c59a674724370321d7d1b914ebc9a46 $
//
#ifndef XCrystalPlanarMoliereElectronDensity_h
#define XCrystalPlanarMoliereElectronDensity_h

#include "XVCrystalPlanarAnalytical.hh"

class XCrystalPlanarMoliereElectronDensity:public XVCrystalPlanarAnalytical {

private:
    G4double fAlfa[3];
    G4double fBeta[3];
  
public:
    G4double ComputeECForSinglePlane(G4double,XPhysicalLattice*);
    
    //Contructors
    XCrystalPlanarMoliereElectronDensity();
    ~XCrystalPlanarMoliereElectronDensity();
};

#endif

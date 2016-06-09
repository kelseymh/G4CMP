/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id: e22ecc58a0c36f7f641db1986900b03d07699a1e $
//
#ifndef XCrystalIntegratedDensityPlanar_h
#define XCrystalIntegratedDensityPlanar_h

#include "XVCrystalIntegratedDensity.hh"

class XCrystalIntegratedDensityPlanar : public XVCrystalIntegratedDensity{

protected:
    virtual G4double ComputeIntegratedDensity(G4double,G4int);

public:   //Contructors
    XCrystalIntegratedDensityPlanar();
    ~XCrystalIntegratedDensityPlanar();
};

#endif

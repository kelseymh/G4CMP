/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id: 6ccd5344828bb56f77d1b6b95b23626281bf18f8 $
//
#ifndef XCrystalCharacteristicArray_h
#define XCrystalCharacteristicArray_h

#include "XVCrystalCharacteristic.hh"

class XCrystalCharacteristicArray: public XVCrystalCharacteristic {

private:
    std::vector<XVCrystalCharacteristic*> fCharacteristicVector;
protected:
    
public:
    std::vector<XVCrystalCharacteristic*> GetCharacteristicVector();
    
    virtual G4ThreeVector ComputeEC(G4ThreeVector,XPhysicalLattice*);
    virtual G4ThreeVector ComputePositionInUnitCell(G4ThreeVector,XPhysicalLattice*);
    
    //Contructors
    XCrystalCharacteristicArray();
    ~XCrystalCharacteristicArray();
};

#endif

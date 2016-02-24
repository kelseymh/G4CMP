/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id: 2d843e7a8ac33144c9794263fd57914ab76ad89b $
//
#ifndef XUnitCell_h
#define XUnitCell_h

#include <iostream>
#include <fstream>
#include <string>
#include "G4ThreeVector.hh"

#include "XLogicalAtomicLattice.hh"
#include "XLogicalAtomicLatticeDiamond.hh"
#include "XLogicalBase.hh"

#define MAXBASENUMBER 32

using namespace std;

class XUnitCell{

private:
    G4int fNumberOfBases;
    XLogicalBase* fBase[MAXBASENUMBER];
    
    G4ThreeVector fSize;
    G4ThreeVector fAngle;
    
public:
    //Retrieval methods
    G4ThreeVector GetSize();
    G4ThreeVector GetAngle();
    XLogicalBase* GetBase(G4int);
    
    
    //Set methods
    void SetSize(G4ThreeVector);
    void SetAngle(G4ThreeVector);
    void SetBase(G4int,XLogicalBase*);
    void AddBase(XLogicalBase*);

    //Calculation methods
    G4double ComputeVolume();

    G4double ComputeMillerOverSizeSquared(G4int,G4int,G4int);
    G4double ComputeMillerPerSizeSquared(G4int,G4int,G4int);

    G4double ComputeReciprocalVectorSquared(G4int,G4int,G4int);
    G4double ComputeReciprocalVector(G4int,G4int,G4int);

    G4double ComputeDirectVectorSquared(G4int,G4int,G4int);
    G4double ComputeDirectVector(G4int,G4int,G4int);
    
    G4double ComputeDirectPeriodSquared(G4int,G4int,G4int);
    G4double ComputeDirectPeriod(G4int,G4int,G4int);

    
    G4double ComputeAtomVolumeDensity();
    G4complex ComputeStructureFactor(G4int,G4int,G4int); //Kittel - chapter 2 Eq. (46)

    //Check method
    G4bool IsOrthogonal();
    G4bool IsCubic();
    
    //Contructors
    XUnitCell();
    ~XUnitCell();
};

#endif

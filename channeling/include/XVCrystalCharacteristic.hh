//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id$
//
#ifndef XVCrystalCharacteristic_h
#define XVCrystalCharacteristic_h

#include "XLatticeManager3.hh"
#include "G4Track.hh"

class XVCrystalCharacteristic {

private:
    XLatticeManager3* fLatticeManager;

protected:
    G4double fThermalVibrationAmplitude; // TO BE MOVED TO XLogicalLattice
    G4ThreeVector fMaximum;
    G4ThreeVector fMinimum;
    
public:
    //retrieval functions
    XPhysicalLattice* GetXPhysicalLattice(G4VPhysicalVolume*);
    XUnitCell* GetXUnitCell(G4VPhysicalVolume*);
    XLogicalLattice* GetLogicalLattice(G4VPhysicalVolume*);

    void SetThermalVibrationAmplitude(G4double);
    G4double GetThermalVibrationAmplitude();
    
    //virtual function to compute value starting from the point in the xtal reference frame and the physical volume of the xtal
    virtual G4ThreeVector ComputeValue(G4ThreeVector,XPhysicalLattice*) = 0;
    virtual G4ThreeVector ComputePositionInUnitCell(G4ThreeVector,XPhysicalLattice*);
    
    virtual G4double ComputeTFScreeningRadius(XPhysicalLattice*);
    
    virtual G4ThreeVector GetMaximum(XPhysicalLattice*);
    virtual G4ThreeVector GetMinimum(XPhysicalLattice*);

    virtual G4ThreeVector ComputeMaximum(XPhysicalLattice*);
    virtual G4ThreeVector ComputeMinimum(XPhysicalLattice*);

    virtual void PrintOnFile(char*,XPhysicalLattice*,G4double = 1) = 0;
    //Contructors
    XVCrystalCharacteristic();
    ~XVCrystalCharacteristic();
};

#endif

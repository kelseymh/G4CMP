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
#ifndef XCrystalIntegratedDensity_h
#define XCrystalIntegratedDensity_h

#include "XVCrystalCharacteristic.hh"

class XCrystalIntegratedDensity {

public:
    void SetIntegrationPoints(unsigned int,unsigned int);
    unsigned int GetIntegrationPoints(unsigned int);
    unsigned int GetIntegrationPoints();

    void SetNumberOfPoints(unsigned int,unsigned int);
    unsigned int GetNumberOfPoints(unsigned int);
    unsigned int GetNumberOfPoints();

    G4double GetStep(unsigned int);

    void SetNucleiDensity(XVCrystalCharacteristic*);
    XVCrystalCharacteristic* GetNucleiDensity();
    
    void SetElectronDensity(XVCrystalCharacteristic*);
    XVCrystalCharacteristic* GetElectronDensity();

    void SetPotential(XVCrystalCharacteristic*);
    XVCrystalCharacteristic* GetPotential();

    void SetXPhysicalLattice(XPhysicalLattice*);
    XPhysicalLattice* GetXPhysicalLattice();
    
    G4double GetValue(G4double);
    
    void InitializeTable();
    G4bool HasBeenInitialized();
    //now it checks only of the table is initialized, it does not check if the particular crystal is initialized. To be changed in the future!
    
    void PrintOnFile(char*);
    
protected:
    G4double ComputeValue(G4double);

private:
    XPhysicalLattice* fLattice;
    XVCrystalCharacteristic* fElectronDensity;
    XVCrystalCharacteristic* fNucleiDensity;
    XVCrystalCharacteristic* fPotential;

    G4double fPotentialMinimum;
    G4double fPotentialRange;
    
    std::vector<G4double> fTable;
    unsigned int fNumberOfPoints[3];
    unsigned int fIntegrationPoints[3];
    

public:   //Contructors
    XCrystalIntegratedDensity(XVCrystalCharacteristic*,XVCrystalCharacteristic*,XVCrystalCharacteristic*);
    ~XCrystalIntegratedDensity();
};

#endif

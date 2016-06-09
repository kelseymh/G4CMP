/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id: 3578c019c94c1f34dd5a3a336f83283629ba702d $
//
#ifndef XCrystalIntegratedDensityHub_h
#define XCrystalIntegratedDensityHub_h

#include "XVCrystalCharacteristic.hh"
#include "XVCrystalIntegratedDensity.hh"
#include "G4ParticleDefinition.hh"

class XCrystalIntegratedDensityHub {

public:
    void SetDensityElectron(XVCrystalCharacteristic*);
    XVCrystalCharacteristic* GetDensityElectron();

    void SetDensityNuclei(XVCrystalCharacteristic*);
    XVCrystalCharacteristic* GetDensityNuclei();

    void SetPotential(XVCrystalCharacteristic*);
    XVCrystalCharacteristic* GetPotential();

    void SetXPhysicalLattice(XPhysicalLattice*);
    XPhysicalLattice* GetXPhysicalLattice();

    void SetIntegratedDensityNuclei(XVCrystalIntegratedDensity*,G4int);
    XVCrystalIntegratedDensity* GetIntegratedDensityNuclei(G4int);

    void SetIntegratedDensityElectron(XVCrystalIntegratedDensity*,G4int);
    XVCrystalIntegratedDensity* GetIntegratedDensityElectron(G4int);

    void PrintOnFiles(const G4String&);
    G4bool HasBeenInitialized(XPhysicalLattice*);

    G4double GetIntegratedDensityElectron(G4double,XPhysicalLattice*,G4int);
    G4double GetIntegratedDensityNuclei(G4double,XPhysicalLattice*,G4int);

public:
    void InitializeTables();
    
private:
    XPhysicalLattice* fLattice;
    XVCrystalCharacteristic* fDensityElectron;
    XVCrystalCharacteristic* fDensityNuclei;
    XVCrystalCharacteristic* fPotential;

    XVCrystalIntegratedDensity* fIntDensElectronPositive;
    XVCrystalIntegratedDensity* fIntDensNucleiPositive;

    XVCrystalIntegratedDensity* fIntDensElectronNegative;
    XVCrystalIntegratedDensity* fIntDensNucleiNegative;

public:   //Contructors
    XCrystalIntegratedDensityHub();
    ~XCrystalIntegratedDensityHub();
};

#endif

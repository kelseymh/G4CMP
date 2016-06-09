/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id: 3cd3b17a2281402405a1118506a79f427d310f51 $
//
#ifndef XVCrystalIntegratedDensity_h
#define XVCrystalIntegratedDensity_h

#include "XVCrystalCharacteristic.hh"
#include "G4PhysicsLinearVector.hh"

class XVCrystalIntegratedDensity {

public:
    void SetIntegrationPoints(unsigned int,unsigned int);
    unsigned int GetIntegrationPoints(unsigned int);
    unsigned int GetIntegrationPoints();

    void SetNumberOfPoints(unsigned int);
    unsigned int GetNumberOfPoints();

    void SetDensity(XVCrystalCharacteristic*);
    XVCrystalCharacteristic* GetDensity();
    
    void SetPotential(XVCrystalCharacteristic*);
    XVCrystalCharacteristic* GetPotential();

    void SetXPhysicalLattice(XPhysicalLattice*);
    XPhysicalLattice* GetXPhysicalLattice();

    void SetParticleCharge(G4int);
    G4int GetParticleCharge();

    void PrintOnFile(const G4String&);
    void ReadFromFile(const G4String&);
    
    G4bool HasBeenInitialized(XPhysicalLattice*,G4int);

    G4double GetIntegratedDensity(G4double,XPhysicalLattice*,G4int);

protected:
    G4double GetStep();
    
    virtual void ComputePotentialParameters();

public:
    virtual void InitializeTable();
    
protected:
    virtual G4double ComputeIntegratedDensity(G4double,G4int);
    G4double FindCatmullRomInterpolate(G4double &p0, G4double &p1, G4double &p2, G4double &p3, G4double &x);

private:
    XPhysicalLattice* fLattice;
    G4int fParticleCharge;
    XVCrystalCharacteristic* fDensity;
    XVCrystalCharacteristic* fPotential;

protected:
    G4double fPotentialMinimum;
    G4double fPotentialMaximum;
    G4double fPotentialRange;

private:
    G4PhysicsLinearVector* fTableVector;
    unsigned int fNumberOfPoints;
    unsigned int fIntegrationPoints[3];
    

public:   //Contructors
    XVCrystalIntegratedDensity();
    ~XVCrystalIntegratedDensity();
};

#endif

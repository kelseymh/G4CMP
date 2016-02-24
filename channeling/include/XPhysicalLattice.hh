/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file exoticphysics/phonon/include/XPhysicalLattice.hh
/// \brief Definition of the XPhysicalLattice class
//
// $Id: 1e69488374636eba6c67ce80a03e52efb9468166 $
//
#ifndef XPhysicalLattice_h
#define XPhysicalLattice_h 1

#include "XLogicalLattice.hh"
#include "G4ThreeVector.hh"
#include "G4AffineTransform.hh"
#include "XUnitCell.hh"

class G4VPhysicalVolume;

class XPhysicalLattice{
    
private:
    G4double fTheta, fPhi, fOmega;
    XLogicalLattice* fLattice;
    G4VPhysicalVolume* fVolume;
    
    double fA;       //Scaling constant for Anh.Dec. mean free path
    double fB;       //Scaling constant for Iso.Scat. mean free path
    double fDosL;    //Density of states for L-phonons
    double fDosST;   //Density of states for ST-phonons
    double fDosFT;    //Density of states for FT-phonons
    double fBeta, fGamma, fLambda, fMu; //dynamical constants for material
    
    
public:
    G4AffineTransform fLocalToGlobal;
    G4AffineTransform fGlobalToLocal;
    
    void SetDynamicalConstants(double, double, double, double);
    void SetScatteringConstant(G4double);
    void SetAnhDecConstant(G4double);
    void SetLDOS(double);
    void SetSTDOS(double);
    void SetFTDOS(double);
    
    double GetBeta();
    double GetGamma();
    double GetLambda();
    double GetMu();
    G4double GetScatteringConstant();
    G4double GetAnhDecConstant();
    double GetLDOS();
    double GetSTDOS();
    double GetFTDOS();
    
    XPhysicalLattice();
    XPhysicalLattice(G4VPhysicalVolume*, XLogicalLattice*);
    ~XPhysicalLattice();
    double MapKtoV(int, G4ThreeVector);
    G4ThreeVector MapKtoVDir(int, G4ThreeVector);
    G4VPhysicalVolume* GetVolume();
    void SetPhysicalVolume(G4VPhysicalVolume*);
    void SetXLogicalLattice(XLogicalLattice*);
    void SetLatticeOrientation(G4double, G4double);
    void SetLatticeOrientation(G4double, G4double, G4double);
    void SetMillerOrientation(int, int, int);
    
    
    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
    // Begin channeling specific code
    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
    
public:
    //set methods
    void SetUnitCell(XUnitCell*);
    G4ThreeVector ProjectMomentumVectorFromWorldToLattice(G4ThreeVector,G4ThreeVector);
    G4ThreeVector ProjectMomentumVectorFromLatticeToWorld(G4ThreeVector,G4ThreeVector);
    G4ThreeVector GetLatticeDirection(G4ThreeVector);

    //retrieval methods
    XUnitCell* GetXUnitCell();
    XLogicalLattice* GetLogicalLattice();
    G4int GetMiller(G4int);
    G4ThreeVector GetLatticeAngles();

    G4ThreeVector GetCurvatureRadius();
    void SetCurvatureRadius(G4ThreeVector);
    G4ThreeVector ComputeBendingAngle(G4ThreeVector);
    
    G4bool IsBent();
    
    //general functions
    G4double ComputeInterplanarPeriod();
    
    void SetThermalVibrationAmplitude(G4double);
    G4double GetThermalVibrationAmplitude();
    void ComputeThermalVibrationAmplitude();

private:
    G4ThreeVector fCurvatureRadius;
    G4double fThermalVibrationAmplitude; // TO BE MOVED TO XLogicalLattice
    G4int fMillerOrientation[3];
    XUnitCell* fUnitCell;
};

#endif

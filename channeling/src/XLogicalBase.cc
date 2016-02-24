/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#include "XLogicalBase.hh"
#include "G4PhysicalConstants.hh"
#include <cmath>

XLogicalBase::XLogicalBase(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XLogicalBase::XLogicalBase(G4Element* vElement,XLogicalAtomicLattice* vLattice){
    SetElement(vElement);
    SetLattice(vLattice);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
XLogicalBase::~XLogicalBase(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XLogicalAtomicLattice* XLogicalBase::GetLattice(){
    return fLattice;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4Element* XLogicalBase::GetElement(){
    return fElement;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XLogicalBase::SetLattice(XLogicalAtomicLattice* vLattice){
    fLattice = vLattice;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XLogicalBase::SetElement(G4Element* vElement){
    fElement = vElement;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XLogicalBase::ComputeAtomicFormFactor(){
    return 1.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4complex XLogicalBase::ComputeStructureFactorSingleAtomicKind(G4int h,G4int k ,G4int l){
    G4double vAtomicFormFactor = ComputeAtomicFormFactor();
    G4complex vResult = GetLattice()->ComputeGeometricalStructureFactorSingleKind(h,k,l);
    vResult = G4complex(vResult.real() * vAtomicFormFactor,vResult.imag() * vAtomicFormFactor);
    return vResult;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



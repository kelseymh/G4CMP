/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#include "XCrystalPlanarMoliereElectronDensity.hh"

XCrystalPlanarMoliereElectronDensity::XCrystalPlanarMoliereElectronDensity(){
    fAlfa[0] = 0.1;
    fAlfa[1] = 0.55;
    fAlfa[2] = 0.35;
    
    fBeta[0] = 6.0;
    fBeta[1] = 1.2;
    fBeta[2] = 0.3;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XCrystalPlanarMoliereElectronDensity::~XCrystalPlanarMoliereElectronDensity(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XCrystalPlanarMoliereElectronDensity::ComputeECForSinglePlane(G4double vXposition,XPhysicalLattice* vLattice){

    G4double aTF = ComputeTFScreeningRadius(vLattice);

    G4double vValueForSinglePlane = 0.;
    for(unsigned int i=0;i<3;i++){
        vValueForSinglePlane += ( fAlfa[i] * fBeta[i] * exp( - fabs(vXposition) * fBeta[i] / aTF ) * ( 1. + fBeta[i] * fabs(vXposition) / aTF) );
    }
    
    vValueForSinglePlane /= aTF;
    
    vValueForSinglePlane *= 0.25;

    return vValueForSinglePlane;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

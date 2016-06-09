/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#include "XCrystalPlanarMoliereElectricField.hh"

XCrystalPlanarMoliereElectricField::XCrystalPlanarMoliereElectricField(){
    fAlfa[0] = 0.1;
    fAlfa[1] = 0.55;
    fAlfa[2] = 0.35;
    
    fBeta[0] = 6.0;
    fBeta[1] = 1.2;
    fBeta[2] = 0.3;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XCrystalPlanarMoliereElectricField::~XCrystalPlanarMoliereElectricField(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XCrystalPlanarMoliereElectricField::ComputeECForSinglePlane(G4double vXposition,XPhysicalLattice* vLattice){

    G4VPhysicalVolume* vVolume = vLattice->GetVolume();

    G4double aTF = ComputeTFScreeningRadius(vLattice);

    G4double vValueForSinglePlane = 0.;
    for(unsigned int i=0;i<3;i++){
        vValueForSinglePlane += ( fAlfa[i] * exp( - fabs(vXposition) * fBeta[i] / aTF ) );
    }

    vValueForSinglePlane *= 2. * M_PI * GetXUnitCell(vVolume)->ComputeDirectPeriod(GetXPhysicalLattice(vVolume)->GetMiller(0),GetXPhysicalLattice(vVolume)->GetMiller(1),GetXPhysicalLattice(vVolume)->GetMiller(2));
    
    vValueForSinglePlane *= (CLHEP::elm_coupling);
    
    vValueForSinglePlane *= (GetXUnitCell(vVolume)->ComputeAtomVolumeDensity());
    
    G4int vSign = +1;
    
    if(vXposition < 0.) vSign = -1;

    vValueForSinglePlane *= vSign;
    
    return vValueForSinglePlane;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

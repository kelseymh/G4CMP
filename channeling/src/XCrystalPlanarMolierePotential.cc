/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#include "XCrystalPlanarMolierePotential.hh"

XCrystalPlanarMolierePotential::XCrystalPlanarMolierePotential(){
    fAlfa[0] = 0.1;
    fAlfa[1] = 0.55;
    fAlfa[2] = 0.35;
    
    fBeta[0] = 6.0;
    fBeta[1] = 1.2;
    fBeta[2] = 0.3;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XCrystalPlanarMolierePotential::~XCrystalPlanarMolierePotential(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XCrystalPlanarMolierePotential::ComputeECForSinglePlane(G4double vXposition,XPhysicalLattice* vLattice){
    
    G4VPhysicalVolume* vVolume = vLattice->GetVolume();
    G4double aTF = ComputeTFScreeningRadius(vLattice);

    G4double vValueForSinglePlane = 0.;
    
    for(unsigned int i=0;i<3;i++){
        vValueForSinglePlane += ( fAlfa[i]/fBeta[i] * exp( - fabs(vXposition) * fBeta[i] / aTF ) );
    }
    
    vValueForSinglePlane *= 2. * M_PI * GetXUnitCell(vVolume)->ComputeDirectPeriod(GetXPhysicalLattice(vVolume)->GetMiller(0),GetXPhysicalLattice(vVolume)->GetMiller(1),GetXPhysicalLattice(vVolume)->GetMiller(2));
    
    vValueForSinglePlane *= aTF;

    vValueForSinglePlane *= (CLHEP::elm_coupling);

    vValueForSinglePlane *= (GetXUnitCell(vVolume)->ComputeAtomVolumeDensity());

    return vValueForSinglePlane;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XCrystalPlanarMolierePotential::ComputeMaximum(XPhysicalLattice* vLattice){

    G4double vMaximum = GetEC(G4ThreeVector(0.,0.,0.),vLattice).y();
    
    return vMaximum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XCrystalPlanarMolierePotential::ComputeMinimum(XPhysicalLattice* vLattice){
    G4VPhysicalVolume* vVolume = vLattice->GetVolume();
    G4double vInterplanarDistance = GetXUnitCell(vVolume)->ComputeDirectPeriod(GetXPhysicalLattice(vVolume)->GetMiller(0),GetXPhysicalLattice(vVolume)->GetMiller(1),GetXPhysicalLattice(vVolume)->GetMiller(2));
    
    G4double vMinimum = GetEC(G4ThreeVector(0.,vInterplanarDistance/2.,0.),vLattice).y();
    
    return vMinimum;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

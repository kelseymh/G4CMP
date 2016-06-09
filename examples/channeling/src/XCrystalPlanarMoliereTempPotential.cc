/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#include "XCrystalPlanarMoliereTempPotential.hh"

XCrystalPlanarMoliereTempPotential::XCrystalPlanarMoliereTempPotential(){
    fAlfa[0] = 0.1;
    fAlfa[1] = 0.55;
    fAlfa[2] = 0.35;
    
    fBeta[0] = 6.0;
    fBeta[1] = 1.2;
    fBeta[2] = 0.3;
    
    for(unsigned int i=0;i<3;i++) {
        fGamma[i] = fAlfa[i]/fBeta[i];
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XCrystalPlanarMoliereTempPotential::~XCrystalPlanarMoliereTempPotential(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XCrystalPlanarMoliereTempPotential::ComputeECForSinglePlane(G4double vX,XPhysicalLattice* vLattice){
    
    G4VPhysicalVolume* vVolume = vLattice->GetVolume();

    G4double aTF = ComputeTFScreeningRadius(vLattice);
    G4double vTVA = vLattice->GetThermalVibrationAmplitude();
    
    G4double vTau[3];
    for(unsigned int i=0;i<3;i++) vTau[i] = (pow( vTVA / aTF * fBeta[i] , 2. ) / 2.0);

    
    G4double vValueForSinglePlane = 0.;
    
    for(unsigned int i=0;i<3;i++){
        G4double vTemp = 0.;
        vTemp += ( exp(-vX/ aTF * fBeta[i] ) * erfc((vTVA / aTF * fBeta[i] - vX/ vTVA) / pow(2.,0.5)) );
        vTemp += ( exp( vX/ aTF * fBeta[i] ) * erfc((vTVA / aTF * fBeta[i] + vX/ vTVA) / pow(2.,0.5)) );
        vValueForSinglePlane += ( vTemp * fGamma[i] * exp( vTau[i] ) /2.0);
    }
    
    vValueForSinglePlane *= 2. * M_PI * GetXUnitCell(vVolume)->ComputeDirectPeriod(GetXPhysicalLattice(vVolume)->GetMiller(0),GetXPhysicalLattice(vVolume)->GetMiller(1),GetXPhysicalLattice(vVolume)->GetMiller(2));
    
    vValueForSinglePlane *= aTF;

    vValueForSinglePlane *= (CLHEP::elm_coupling);

    vValueForSinglePlane *= (GetXUnitCell(vVolume)->ComputeAtomVolumeDensity());

    return vValueForSinglePlane;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XCrystalPlanarMoliereTempPotential::ComputeMaximum(XPhysicalLattice* vLattice){

    G4double vMaximum = GetEC(G4ThreeVector(0.,0.,0.),vLattice).x();
    
    return vMaximum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XCrystalPlanarMoliereTempPotential::ComputeMinimum(XPhysicalLattice* vLattice){
    G4VPhysicalVolume* vVolume = vLattice->GetVolume();
    G4double vInterplanarDistance = GetXUnitCell(vVolume)->ComputeDirectPeriod(GetXPhysicalLattice(vVolume)->GetMiller(0),GetXPhysicalLattice(vVolume)->GetMiller(1),GetXPhysicalLattice(vVolume)->GetMiller(2));
    
    G4double vMinimum = GetEC(G4ThreeVector(vInterplanarDistance/2.,0.,0.),vLattice).x();
    
    return vMinimum;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#include "XCrystalPlanarNucleiDensity.hh"

XCrystalPlanarNucleiDensity::XCrystalPlanarNucleiDensity(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XCrystalPlanarNucleiDensity::~XCrystalPlanarNucleiDensity(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XCrystalPlanarNucleiDensity::ComputeECForSinglePlane(G4double vXposition,XPhysicalLattice* vLattice){
    
    G4double vValueForSinglePlane = exp( - 0.5 * pow(vXposition/vLattice->GetThermalVibrationAmplitude(),2.0 ) );

    vValueForSinglePlane /= (vLattice->GetThermalVibrationAmplitude());
    vValueForSinglePlane /= ( sqrt( 2 * M_PI) );
    
    return vValueForSinglePlane;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

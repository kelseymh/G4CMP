/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#include "XCrystalIntegratedDensityPlanar.hh"

XCrystalIntegratedDensityPlanar::XCrystalIntegratedDensityPlanar(){
    SetNumberOfPoints(256);
    
    SetIntegrationPoints(0,1024);
    SetIntegrationPoints(1,1);
    SetIntegrationPoints(2,1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XCrystalIntegratedDensityPlanar::~XCrystalIntegratedDensityPlanar(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XCrystalIntegratedDensityPlanar::ComputeIntegratedDensity(G4double vPotentialInitial, G4int vParticleCharge){
    
    unsigned int i1 = 0;
    
    G4ThreeVector vPositionTemp = G4ThreeVector(0.,0.,0.);
    G4double vInterplanarPeriod = GetXPhysicalLattice()->ComputeInterplanarPeriod();
    G4double vPotential = 0.;
    G4double vDensity = 0.;
    G4double xPos = 0.;
    G4double xMin = +vInterplanarPeriod;
    G4double xMax = -vInterplanarPeriod;
    
    while(i1<GetIntegrationPoints(0)){
        
        xPos = G4double(G4double(i1)/G4double(GetIntegrationPoints(0))*vInterplanarPeriod);
        
        vPositionTemp.setX(xPos + vInterplanarPeriod / 2.);
        
        vPotential = G4double(vParticleCharge) * GetPotential()->GetEC(vPositionTemp,GetXPhysicalLattice()).x();
        
        if(vPotential <= vPotentialInitial){
            vDensity += GetDensity()->GetEC(vPositionTemp,GetXPhysicalLattice()).x();
            if(xPos < xMin){
                xMin = xPos;
            }
            if(xPos > xMax){
                xMax = xPos;
            }
        }
        
        i1++;
    };
    
    
    vDensity *= vInterplanarPeriod; // time correction for the density

    vDensity *= vInterplanarPeriod/fabs(xMax - xMin); // time correction for the density

    vDensity /= GetIntegrationPoints(0);
    
    return vDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

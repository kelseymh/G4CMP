//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//

#include "XVCrystalIntegratedDensity.hh"

XVCrystalIntegratedDensity::XVCrystalIntegratedDensity(){
    fTable.clear();
    
    fNumberOfPoints = 512;
    fIntegrationPoints[0] = 32;
    fIntegrationPoints[1] = 32;
    fIntegrationPoints[2] = 32;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XVCrystalIntegratedDensity::~XVCrystalIntegratedDensity(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XVCrystalIntegratedDensity::SetIntegrationPoints(unsigned int vIndex,unsigned int vIntegrationPoints){
    if(vIndex<3) {
        if(vIntegrationPoints > 0){
            fIntegrationPoints[vIndex] = vIntegrationPoints;
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

unsigned int XVCrystalIntegratedDensity::GetIntegrationPoints(unsigned int vIndex){
    if(vIndex<3) {
        return fIntegrationPoints[vIndex];
    }
    else{
        return 0;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

unsigned int XVCrystalIntegratedDensity::GetIntegrationPoints(){
    return fIntegrationPoints[0]*fIntegrationPoints[1]*fIntegrationPoints[2];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XVCrystalIntegratedDensity::SetNumberOfPoints(unsigned int vNumberOfPoints){
    fNumberOfPoints = vNumberOfPoints;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

unsigned int XVCrystalIntegratedDensity::GetNumberOfPoints(){
    return fNumberOfPoints;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XVCrystalIntegratedDensity::SetPotential(XVCrystalCharacteristic* vPotential){
    fPotential = vPotential;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XVCrystalCharacteristic* XVCrystalIntegratedDensity::GetPotential(){
    return fPotential;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XVCrystalIntegratedDensity::SetDensity(XVCrystalCharacteristic* vDensity){
    fDensity = vDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XVCrystalCharacteristic* XVCrystalIntegratedDensity::GetDensity(){
    return fDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XVCrystalIntegratedDensity::SetXPhysicalLattice(XPhysicalLattice* vLattice){
    fLattice = vLattice;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XPhysicalLattice*  XVCrystalIntegratedDensity::GetXPhysicalLattice(){
    return fLattice;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XVCrystalIntegratedDensity::GetStep(){
    return fPotentialRange / fNumberOfPoints;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool XVCrystalIntegratedDensity::HasBeenInitialized(XPhysicalLattice* vLattice){
    if(fTable.size()==0) return false;
    if(vLattice!=GetXPhysicalLattice()) return false;
    else return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XVCrystalIntegratedDensity::InitializeTable(){
    fPotentialMinimum = fPotential->GetMinimum(fLattice);
    fPotentialRange = fPotential->GetMaximum(fLattice) - fPotentialMinimum;
    
    G4cout << "XVCrystalIntegratedDensity::InitializeTable()::Potential Range =  " << fPotentialRange/eV << std::endl;
    
    for(unsigned int i=0;i<GetNumberOfPoints();i++){
        fTable.push_back(ComputeValue(fPotentialMinimum + fPotentialRange * G4double(i+1) / G4double(fNumberOfPoints)));
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XVCrystalIntegratedDensity::GetValue(G4double vPotential,XPhysicalLattice* vLattice){
    if(!HasBeenInitialized(vLattice)) return -1.;
    else if(vPotential >= (fPotentialMinimum+fPotentialRange)) return 1.;
    else if(vPotential <= fPotentialMinimum) return 0.;
    else{
        double vX = (vPotential - fPotentialMinimum) / fPotentialRange ;
        unsigned int vIndex = round(vX * double(fNumberOfPoints));
        
        bool bIndexNearLeft = false;
        if(vIndex == trunc(vX * double(fNumberOfPoints))) bIndexNearLeft = true;
                
        if( (vIndex > 2) && (vIndex < (fNumberOfPoints - 2)) ){
            if(bIndexNearLeft){
                return FindCatmullRomInterpolate(fTable[vIndex-1],fTable[vIndex],fTable[vIndex+1],fTable[vIndex+2],vX);
            }
            else{
                return FindCatmullRomInterpolate(fTable[vIndex-2],fTable[vIndex-1],fTable[vIndex],fTable[vIndex+1],vX);
            }
        }
        else{
            if(vIndex <= 2){
                return FindCatmullRomInterpolate(fTable[0],fTable[1],fTable[2],fTable[3],vX);
            }
            else if(vIndex >=  (fNumberOfPoints - 2)){
                return FindCatmullRomInterpolate(fTable[fNumberOfPoints-3],fTable[fNumberOfPoints-2],fTable[fNumberOfPoints-1],fTable[fNumberOfPoints],vX);
            }
            else{
                return 6.;

                return fTable[vIndex];
            }
        }
    }
    return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XVCrystalIntegratedDensity::ComputeValue(G4double vPotentialInitial){
    
    unsigned int i1,i2,i3;
    i1 = i2 = i3 = 0;
    
    G4ThreeVector vPositionTemp = G4ThreeVector(0.,0.,0.);
    G4double vDensity = 0.;
    
    G4ThreeVector vSize = fLattice->GetXUnitCell()->GetSize();
    while(i1<fIntegrationPoints[2]){
        vPositionTemp.setY(G4double(G4double(i3)/G4double(fIntegrationPoints[2])*vSize.z()));
        while(i1<fIntegrationPoints[1]){
            vPositionTemp.setZ(G4double(G4double(i2)/G4double(fIntegrationPoints[1])*vSize.y()));
            while(i1<fIntegrationPoints[0]){
                vPositionTemp.setX(G4double(G4double(i1)/G4double(fIntegrationPoints[0])*vSize.x()));
                if(fPotential->ComputeValue(vPositionTemp,fLattice).x() < vPotentialInitial){
                    vDensity += fDensity->ComputeValue(vPositionTemp,fLattice).x();
                }
                i1++;
            };
            i2++;
        };
        i3++;
    };
    
    vDensity *= fLattice->GetXUnitCell()->ComputeVolume();
    vDensity /= GetIntegrationPoints();
    
    return vDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XVCrystalIntegratedDensity::FindCatmullRomInterpolate(G4double &p0,G4double &p1, G4double &p2,G4double &p3,G4double &x)
{
    double a0, a1 , a2 , a3 , x2;
    
    x2 = x * x;
    a0 = -0.5 * p0 + 1.5 * p1 - 1.5 * p2 + 0.5 * p3;
    a1 = p0 - 2.5 * p1 + 2.0 * p2 - 0.5 * p3;
    a2 = -0.5 * p0 + 0.5 * p2;
    a3 = p1;
    
    return (a0 * x * x2 + a1 * x2 + a2 * x + a3);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XVCrystalIntegratedDensity::PrintOnFile(const G4String& filename,XPhysicalLattice* vLattice){
    std::ofstream vFileOut;
    vFileOut.open(filename);
    G4double vStep = GetStep() / eV;
    
    vFileOut << "energy,dens" << std::endl;
    for(unsigned int i = 0;i < fNumberOfPoints;i++){
        vFileOut << i * vStep << "," << fTable[i] << std::endl;
    }
    vFileOut.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

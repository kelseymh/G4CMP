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

#include "XThomasFermiScreeningRadius.hh"

XThomasFermiScreeningRadius::XThomasFermiScreeningRadius(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XThomasFermiScreeningRadius::~XThomasFermiScreeningRadius(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XThomasFermiScreeningRadius::ComputeTFScreeningRadius(const G4Track& aTrack){
        
    G4double vTFSR = Bohr_radius * 0.88534;
    
    if(aTrack.GetParticleDefinition()->GetParticleName() == "proton"){
        vTFSR /= (std::pow(aTrack.GetMaterial()->GetZ(),0.333333333));
    }
    else{
        vTFSR /= (std::pow(aTrack.GetMaterial()->GetZ(),0.23) + std::pow(aTrack.GetParticleDefinition()->GetPDGCharge(),0.23));
    }
    
    return vTFSR;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
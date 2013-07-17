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
// $Id: DetectorConstruction.hh,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

#include "ChannelingParticleUserInfo.hh"

ChannelingParticleUserInfo::ChannelingParticleUserInfo(){
    bHasBeenUnderCoherentEffect = 0;
    
    fNucleiDensity = 1.0;
    fNucleiDensityPreviousStep = 1.0;
    
    fElectronDensity = 1.0;
    fElectronDensityPreviousStep = 1.0;
    
    fNumberOfDechanneling = 0;
    
    fMomentumInChanneling = G4ThreeVector(0.,0.,0.);
    fPositionInChanneling = G4ThreeVector(0.,0.,0.);
    
    fMomentumInChannelingInitial = G4ThreeVector(DBL_MAX,DBL_MAX,DBL_MAX);
    fPositionInChannelingInitial = G4ThreeVector(DBL_MAX,DBL_MAX,DBL_MAX);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ChannelingParticleUserInfo::~ChannelingParticleUserInfo(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChannelingParticleUserInfo::SetCoherentEffect(G4int flag){
    bHasBeenUnderCoherentEffect = flag;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int ChannelingParticleUserInfo::HasBeenUnderCoherentEffect(){
    return bHasBeenUnderCoherentEffect;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChannelingParticleUserInfo::SetNucleiDensity(G4double density){
    fNucleiDensity = density;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double ChannelingParticleUserInfo::GetNucleiDensity(){
    return fNucleiDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChannelingParticleUserInfo::SetElectronDensity(G4double density){
    fElectronDensity = density;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double ChannelingParticleUserInfo::GetElectronDensity(){
    return fElectronDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double ChannelingParticleUserInfo::GetNucleiDensityPreviousStep(){
    return fNucleiDensityPreviousStep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double ChannelingParticleUserInfo::GetElectronDensityPreviousStep(){
    return fElectronDensityPreviousStep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChannelingParticleUserInfo::StoreDensityPreviousStep(){
    fElectronDensityPreviousStep = fElectronDensity;
    fNucleiDensityPreviousStep = fNucleiDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector ChannelingParticleUserInfo::GetMomentumChanneled(){
    return fMomentumInChanneling;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChannelingParticleUserInfo::SetMomentumChanneled(G4ThreeVector momentum){
    fMomentumInChanneling = momentum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector ChannelingParticleUserInfo::GetPositionChanneled(){
    return fPositionInChanneling;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChannelingParticleUserInfo::SetPositionChanneled(G4ThreeVector position){
    fPositionInChanneling = position;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector ChannelingParticleUserInfo::GetMomentumChanneledInitial(){
    return fMomentumInChannelingInitial;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChannelingParticleUserInfo::SetMomentumChanneledInitial(G4ThreeVector momentum){
    fMomentumInChannelingInitial = momentum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector ChannelingParticleUserInfo::GetPositionChanneledInitial(){
    return fPositionInChannelingInitial;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChannelingParticleUserInfo::SetEnergyChanneled(G4double energy){
    fEnergyInChanneling = energy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double ChannelingParticleUserInfo::GetEnergyChanneled(){
    return fEnergyInChanneling;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChannelingParticleUserInfo::SetEnergyChanneledInitial(G4double energy){
    fEnergyInChannelingInitial = energy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double ChannelingParticleUserInfo::GetEnergyChanneledInitial(){
    return fEnergyInChannelingInitial;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChannelingParticleUserInfo::SetPositionChanneledInitial(G4ThreeVector position){
    fPositionInChannelingInitial = position;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int ChannelingParticleUserInfo::GetNumberOfDechanneling(){
    return fNumberOfDechanneling;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChannelingParticleUserInfo::IncreaseNumberOfDechanneling(){
    fNumberOfDechanneling++;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

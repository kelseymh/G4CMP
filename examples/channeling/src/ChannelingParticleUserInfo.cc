/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id: f695a97eba05391d7248f7d3e0c94819e25529d9 $

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

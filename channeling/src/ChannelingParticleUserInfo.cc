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

ChannelingParticleUserInfo::ChannelingParticleUserInfo()
{
    channelingFlag = false;
    channelingFactor = 1.0;
    preStepChannelingFactor = 1.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ChannelingParticleUserInfo::~ChannelingParticleUserInfo()
{
    channelingFlag = false;
    channelingFactor = 1.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChannelingParticleUserInfo::SetChanneling(bool flag)
{
    channelingFlag = flag;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool ChannelingParticleUserInfo::GetChanneling()
{
    return channelingFlag;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChannelingParticleUserInfo::SetChannelingFactor(G4double factor)
{
    preStepChannelingFactor = channelingFactor;
    channelingFactor = factor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double ChannelingParticleUserInfo::GetChannelingFactor()
{
    return channelingFactor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double ChannelingParticleUserInfo::GetPreStepChannelingFactor()
{
    return preStepChannelingFactor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector ChannelingParticleUserInfo::GetMomentumChanneled(){
    return momentumChanneled;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChannelingParticleUserInfo::SetMomentumChanneled(G4ThreeVector momentum){
    momentumChanneled = momentum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector ChannelingParticleUserInfo::GetPositionChanneled(){
    return positionChanneled;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChannelingParticleUserInfo::SetPositionChanneled(G4ThreeVector position){
    positionChanneled = position;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

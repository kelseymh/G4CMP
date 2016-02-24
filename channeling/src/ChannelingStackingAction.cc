/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file hadronic/Hadr01/src/ChannelingStackingAction.cc
/// \brief Implementation of the ChannelingStackingAction class
//
// $Id: 5588afdfb57082956341332feed61892dc8b5414 $
//
//

#include "ChannelingStackingAction.hh"
#include "ChannelingStackingMessenger.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ChannelingStackingAction::ChannelingStackingAction()
{
    fStackMessenger = new ChannelingStackingMessenger(this);
    fKillSecondary  = false;
    fParticle       = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ChannelingStackingAction::~ChannelingStackingAction()
{
    delete fStackMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
ChannelingStackingAction::ClassifyNewTrack(const G4Track* aTrack)
{
    G4ClassificationOfNewTrack status = fUrgent;
    
    const G4ParticleDefinition* part = aTrack->GetDefinition();
    
    //stack or delete secondaries
    if(aTrack->GetTrackID() > 1) {
        if (fKillSecondary || fParticle == part){
            status = fKill;
        }
        if (part->GetPDGCharge()==0){
            status = fKill;
        }
    }
    return status;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChannelingStackingAction::SetKillStatus(G4bool value)
{
    fKillSecondary = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChannelingStackingAction::SetKill(const G4String& name)
{
    fParticle = G4ParticleTable::GetParticleTable()->FindParticle(name);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

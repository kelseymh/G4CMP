/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file analysis/A01/src/ChannelingDriftChamber.cc
/// \brief Implementation of the ChannelingDriftChamber class
//
// $Id: 7d9f3e2326056890bd85415476ca9d890187330d $
// --------------------------------------------------------------
//
#include "ChannelingDriftChamber.hh"
#include "ChannelingDriftHit.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4Navigator.hh"
#include "G4ios.hh"

ChannelingDriftChamber::ChannelingDriftChamber(G4String name):G4VSensitiveDetector(name)
{
    G4String HCname;
    collectionName.insert(HCname="collection");
    fHCID = -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ChannelingDriftChamber::~ChannelingDriftChamber(){
    ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ChannelingDriftChamber::Initialize(G4HCofThisEvent*HCE)
{
    fHitsCollection = new ChannelingDriftHitsCollection(SensitiveDetectorName,collectionName[0]);
    if(fHCID<0){
        fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
    }
    HCE->AddHitsCollection(fHCID,fHitsCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool ChannelingDriftChamber::ProcessHits(G4Step*aStep,G4TouchableHistory* /*ROhist*/)
{
    G4double charge = aStep->GetTrack()->GetDefinition()->GetPDGCharge();
//    if(charge!=-1) return true;
//    if(aStep->GetTrack()->GetDefinition()->GetParticleName()!="e-") return true;
    if(charge==0) return true;
  
    G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
    G4StepPoint* postStepPoint = aStep->GetPostStepPoint();

    if(!(postStepPoint->GetStepStatus() == fGeomBoundary)) return true;
    
    G4TouchableHistory* theTouchable = (G4TouchableHistory*)(preStepPoint->GetTouchable());
    G4VPhysicalVolume* thePhysical = theTouchable->GetVolume(0); // mother
    G4int copyNo = thePhysical->GetCopyNo();
    G4ThreeVector worldPos = preStepPoint->GetPosition();
    G4ThreeVector localPos = theTouchable->GetHistory()->GetTopTransform().TransformPoint(worldPos);
    
    ChannelingDriftHit* aHit = new ChannelingDriftHit(copyNo);
    aHit->SetLayerID(copyNo);
    aHit->SetWorldPos(worldPos);
    aHit->SetLocalPos(localPos);
    aHit->SetTime(preStepPoint->GetGlobalTime());
    aHit->SetEnergy(preStepPoint->GetTotalEnergy());
    
    fHitsCollection->insert(aHit);
    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ChannelingDriftChamber::EndOfEvent(G4HCofThisEvent* /*HCE*/)
{
    ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#include "G4CMPElectrodeSensitivity.hh"
#include "G4CMPElectrodeHit.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4Navigator.hh"
#include "G4ios.hh"

G4CMPElectrodeHitsCollection* G4CMPElectrodeSensitivity::hitsCollection = NULL;

G4CMPElectrodeSensitivity::G4CMPElectrodeSensitivity(G4String name)
:G4VSensitiveDetector(name) {
  G4String HCname;
  collectionName.insert(HCname="G4CMPElectrodeHit");
  HCID = -1;
}

G4CMPElectrodeSensitivity::~G4CMPElectrodeSensitivity() {}

G4CMPElectrodeHitsCollection* G4CMPElectrodeSensitivity::getHitsCollection() {
  return hitsCollection;
}

void G4CMPElectrodeSensitivity::Initialize(G4HCofThisEvent*HCE) {
  hitsCollection = new G4CMPElectrodeHitsCollection
                   (SensitiveDetectorName,collectionName[0]);
  if(HCID<0)
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(hitsCollection);
  HCE->AddHitsCollection(HCID,hitsCollection);
}

G4bool G4CMPElectrodeSensitivity::ProcessHits(G4Step* aStep,
                                              G4TouchableHistory* /*ROhist*/) {
  G4StepPoint* postStepPoint = aStep->GetPostStepPoint();
  if(postStepPoint->GetStepStatus()==fGeomBoundary &&
                    aStep->GetNonIonizingEnergyDeposit()) {
    G4Track* track = aStep->GetTrack();
    G4int trackID = track->GetTrackID();
    G4String name = track->GetDefinition()->GetParticleName();
    G4double startE = track->GetVertexKineticEnergy();
    G4double startTime = track->GetGlobalTime()-track->GetLocalTime();
    G4double finalTime = track->GetGlobalTime();
    G4double edp = aStep->GetNonIonizingEnergyDeposit();

    G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
    G4VPhysicalVolume* pVol = preStepPoint->GetPhysicalVolume();
    G4AffineTransform toLocal = G4AffineTransform(pVol->GetRotation(),
                                                  pVol->GetTranslation()).Inverse();

    G4ThreeVector startPosition = toLocal.TransformPoint(track->GetVertexPosition());
    G4ThreeVector finalPosition = toLocal.TransformPoint(postStepPoint->GetPosition());

    G4CMPElectrodeHit* aHit = new G4CMPElectrodeHit();
    aHit->SetTrackID(trackID);
    aHit->SetParticleName(name);
    aHit->SetStartTime(startTime);
    aHit->SetFinalTime(finalTime);
    aHit->SetStartEnergy(startE);
    aHit->SetEnergyDeposit(edp);
    aHit->SetStartPosition(startPosition);
    aHit->SetFinalPosition(finalPosition);

    hitsCollection->insert(aHit);
  }
  
  return true;
}

void G4CMPElectrodeSensitivity::EndOfEvent(G4HCofThisEvent* /*HCE*/) {}


#include "AlminumElectrodeSensitivity.hh"

#include "AlminumElectrodeHit.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4Navigator.hh"
#include "G4ios.hh"
#include "Phonon.hh"

using namespace std;

AlminumElectrodeHitsCollection* AlminumElectrodeSensitivity::hitsCollection = NULL;

AlminumElectrodeSensitivity::AlminumElectrodeSensitivity(G4String name)
:G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="AlminumElectrodeHit");
  HCID = -1;
  writer.open("caustic.ssv",fstream::in | fstream::out | fstream::ate);
  writer2.open("timing.ssv", fstream::in | fstream::out | fstream::ate);
}

AlminumElectrodeSensitivity::~AlminumElectrodeSensitivity(){
  writer.close();
  writer2.close();
}

AlminumElectrodeHitsCollection* AlminumElectrodeSensitivity::getHitsCollection(){
  return hitsCollection;
}

void AlminumElectrodeSensitivity::Initialize(G4HCofThisEvent*HCE)
{
  hitsCollection = new AlminumElectrodeHitsCollection
                   (SensitiveDetectorName,collectionName[0]);
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(hitsCollection); }
  HCE->AddHitsCollection(HCID,hitsCollection);
}

G4bool AlminumElectrodeSensitivity::ProcessHits(G4Step* aStep,G4TouchableHistory* /*ROhist*/)
{
  //if(aStep->GetTrack()->GetDefinition()!=Phonon::PhononDefinition()) return true;
  G4double edp = aStep->GetNonIonizingEnergyDeposit();
  //G4double edp = aStep->GetTotalEnergyDeposit();
  if(edp==0.) return true;
  //G4cout << "Energy of Hit = " << aStep->GetTotalEnergyDeposit() << G4endl;
  //G4cout << "NonIonizingEnergy of Hit = " << edp/eV << G4endl;

  G4double charge = aStep->GetTrack()->GetParticleDefinition()->GetPDGCharge(); //Added by Rob for FETSim 9/27/12
  G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
  G4StepPoint* postStepPoint = aStep->GetPostStepPoint();
  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(preStepPoint->GetTouchable());
  G4ThreeVector worldPos = postStepPoint->GetPosition();
  G4ThreeVector localPos
    = theTouchable->GetHistory()->GetTopTransform().TransformPoint(worldPos);

  AlminumElectrodeHit* aHit = new AlminumElectrodeHit();
  aHit->SetTime(postStepPoint->GetGlobalTime());
  aHit->SetEDep(edp);
  aHit->SetWorldPos(worldPos);
  aHit->SetLocalPos(localPos);
  aHit->SetCharge(charge); //Added by Rob for FETSim 9/27/12

  hitsCollection->insert(aHit);
  
  writer<<"\n"<<worldPos.getX()/mm<<","<<worldPos.getY()/mm;
  writer2<<"\n"<<postStepPoint->GetGlobalTime()/ns<<" "<<aHit->GetEDep()/eV;

  return true;
}

void AlminumElectrodeSensitivity::EndOfEvent(G4HCofThisEvent* /*HCE*/)
{;}


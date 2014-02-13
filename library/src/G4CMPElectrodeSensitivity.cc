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

using namespace std;

G4CMPElectrodeHitsCollection* G4CMPElectrodeSensitivity::hitsCollection = NULL;

G4CMPElectrodeSensitivity::G4CMPElectrodeSensitivity(G4String name)
:G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="G4CMPElectrodeHit");
  HCID = -1;
  writer.open("caustic.ssv",fstream::in | fstream::out | fstream::ate);
  writer2.open("timing.ssv", fstream::in | fstream::out | fstream::ate);
}

G4CMPElectrodeSensitivity::~G4CMPElectrodeSensitivity(){
  writer.close();
  writer2.close();
}

G4CMPElectrodeHitsCollection* G4CMPElectrodeSensitivity::getHitsCollection(){
  return hitsCollection;
}

void G4CMPElectrodeSensitivity::Initialize(G4HCofThisEvent*HCE)
{
  hitsCollection = new G4CMPElectrodeHitsCollection
                   (SensitiveDetectorName,collectionName[0]);
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(hitsCollection); }
  HCE->AddHitsCollection(HCID,hitsCollection);
}

G4bool G4CMPElectrodeSensitivity::ProcessHits(G4Step* aStep,G4TouchableHistory* /*ROhist*/)
{
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

  G4CMPElectrodeHit* aHit = new G4CMPElectrodeHit();
  aHit->SetTime(postStepPoint->GetGlobalTime());
  aHit->SetEDep(edp);
  aHit->SetWorldPos(worldPos);
  aHit->SetLocalPos(localPos);
  aHit->SetCharge(charge); //Added by Rob for FETSim 9/27/12

  hitsCollection->insert(aHit);
  
  writer<<"\n"<<worldPos.getX()/mm<<","<<worldPos.getY()/mm;
  writer2<<"\n"<<postStepPoint->GetGlobalTime()/ns<<" "<<aHit->GetEDep()/eV << " " << postStepPoint->GetVelocity()/m*s;
  //writer2<<"\n"<<postStepPoint->GetGlobalTime()/ns<<" "<< aStep->GetPostStepPoint()->GetKineticEnergy()/eV;

  return true;
}

void G4CMPElectrodeSensitivity::EndOfEvent(G4HCofThisEvent* /*HCE*/)
{;}


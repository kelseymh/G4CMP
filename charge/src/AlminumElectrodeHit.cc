#include "AlminumElectrodeHit.hh"

#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

#include "G4ios.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UnitsTable.hh"
#include "G4VisAttributes.hh"

G4Allocator<AlminumElectrodeHit> AlminumElectrodeHitAllocator;

AlminumElectrodeHit::AlminumElectrodeHit()
{
  time = 0.;
  edep = 0.;
}

AlminumElectrodeHit::~AlminumElectrodeHit()
{;}

AlminumElectrodeHit::AlminumElectrodeHit(const AlminumElectrodeHit &right)
: G4VHit() {
  time = right.time;
  edep = right.edep;
  worldPos = right.worldPos;
  localPos = right.localPos;
}

const AlminumElectrodeHit& AlminumElectrodeHit::operator=(const AlminumElectrodeHit &right)
{
  time = right.time;
  edep = right.edep;
  worldPos = right.worldPos;
  localPos = right.localPos;
  return *this;
}

int AlminumElectrodeHit::operator==(const AlminumElectrodeHit &/*right*/) const
{
  return 0;
}

void AlminumElectrodeHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(worldPos);
    circle.SetScreenSize(15);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(0.65,0.65,0.);
    G4VisAttributes attribs(colour);
    attribs.SetStartTime(time);
    attribs.SetEndTime(time+1*ms);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

const std::map<G4String,G4AttDef>* AlminumElectrodeHit::GetAttDefs() const
{
  G4bool isNew;
  std::map<G4String,G4AttDef>* store
    = G4AttDefStore::GetInstance("AlminumElectrodeHit",isNew);
  if (isNew) {
    G4String HitType("HitType");
    (*store)[HitType] = G4AttDef(HitType,"Hit Type","Physics","","G4String");

    G4String Time("Time");
    (*store)[Time] = G4AttDef(Time,"Time","Physics","G4BestUnit","G4double");

    G4String EDep("EDep");
    (*store)[EDep] = G4AttDef(Time,"EDep","Physics","G4BestUnit","G4double");

    G4String Pos("Pos");
    (*store)[Pos] = G4AttDef(Pos, "Position",
		      "Physics","G4BestUnit","G4ThreeVector");
  }
  return store;
}

std::vector<G4AttValue>* AlminumElectrodeHit::CreateAttValues() const
{
  std::vector<G4AttValue>* values = new std::vector<G4AttValue>;

  values->push_back(G4AttValue("HitType","AlminumElectrodeHit",""));

  values->push_back
    (G4AttValue("Time",G4BestUnit(time,"Time"),""));

  values->push_back
    (G4AttValue("EDep",G4BestUnit(edep,"Energy"),""));

  values->push_back
    (G4AttValue("Pos",G4BestUnit(worldPos,"Length"),""));

  return values;
}

void AlminumElectrodeHit::Print()
{
  G4cout << "  time " << time/ns << " (nsec) : at " << localPos
         << "  -- edep = " << edep/eV << " [eV]" << G4endl;
}



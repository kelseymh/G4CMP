#include "G4CMPElectrodeHit.hh"

#include "G4ios.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UnitsTable.hh"
#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"


G4Allocator<G4CMPElectrodeHit> G4CMPElectrodeHitAllocator;

G4CMPElectrodeHit::G4CMPElectrodeHit()
{
  time = 0.;
  edep = 0.;
}

G4CMPElectrodeHit::~G4CMPElectrodeHit()
{;}

G4CMPElectrodeHit::G4CMPElectrodeHit(const G4CMPElectrodeHit &right)
: G4VHit() {
  time = right.time;
  edep = right.edep;
  worldPos = right.worldPos;
  localPos = right.localPos;
}

const G4CMPElectrodeHit& G4CMPElectrodeHit::operator=(const G4CMPElectrodeHit &right)
{
  time = right.time;
  edep = right.edep;
  worldPos = right.worldPos;
  localPos = right.localPos;
  return *this;
}

int G4CMPElectrodeHit::operator==(const G4CMPElectrodeHit &/*right*/) const
{
  return 0;
}

void G4CMPElectrodeHit::Draw()
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

const std::map<G4String,G4AttDef>* G4CMPElectrodeHit::GetAttDefs() const
{
  G4bool isNew;
  std::map<G4String,G4AttDef>* store
    = G4AttDefStore::GetInstance("G4CMPElectrodeHit",isNew);
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

std::vector<G4AttValue>* G4CMPElectrodeHit::CreateAttValues() const
{
  std::vector<G4AttValue>* values = new std::vector<G4AttValue>;

  values->push_back(G4AttValue("HitType","G4CMPElectrodeHit",""));

  values->push_back
    (G4AttValue("Time",G4BestUnit(time,"Time"),""));

  values->push_back
    (G4AttValue("EDep",G4BestUnit(edep,"Energy"),""));

  values->push_back
    (G4AttValue("Pos",G4BestUnit(worldPos,"Length"),""));

  return values;
}

void G4CMPElectrodeHit::Print()
{
  G4cout << "  time " << time/ns << " (nsec) : at " << localPos
         << "  -- edep = " << edep/eV << " [eV]" << G4endl;
}



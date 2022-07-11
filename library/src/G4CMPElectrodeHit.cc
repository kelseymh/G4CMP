/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// 20220222  G4CMP-289 -- Thread-local allocator must be a pointer.

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


G4ThreadLocal G4Allocator<G4CMPElectrodeHit>* G4CMPElectrodeHitAllocator=0;

G4CMPElectrodeHit::G4CMPElectrodeHit() {}

int G4CMPElectrodeHit::operator==(const G4CMPElectrodeHit &/*right*/) const {
  return 0;
}

void G4CMPElectrodeHit::Draw() {
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(finalPos);
    circle.SetScreenSize(15);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(0.65,0.65,0.);
    G4VisAttributes attribs(colour);
    attribs.SetStartTime(finalTime);
    attribs.SetEndTime(finalTime+1*ms);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

const std::map<G4String,G4AttDef>* G4CMPElectrodeHit::GetAttDefs() const {
  G4bool isNew;
  std::map<G4String,G4AttDef>* store
    = G4AttDefStore::GetInstance("G4CMPElectrodeHit",isNew);
  if (isNew) {
    G4String HitType("HitType");
    (*store)[HitType] = G4AttDef(HitType,"Hit Type","Physics","","G4String");

    G4String Time("Time");
    (*store)[Time] = G4AttDef(Time,"Time","Physics","G4BestUnit","G4double");

    G4String EDepName("EDep");
    (*store)[EDepName] = G4AttDef(EDepName, "EDep","Physics","G4BestUnit","G4double");

    G4String Pos("Pos");
    (*store)[Pos] = G4AttDef(Pos, "Position",
		      "Physics","G4BestUnit","G4ThreeVector");
  }
  return store;
}

std::vector<G4AttValue>* G4CMPElectrodeHit::CreateAttValues() const {
  std::vector<G4AttValue>* values = new std::vector<G4AttValue>;

  values->push_back(G4AttValue("HitType","G4CMPElectrodeHit",""));

  values->push_back
    (G4AttValue("Time",G4BestUnit(finalTime,"Time"),""));

  values->push_back
    (G4AttValue("EDep",G4BestUnit(EDep,"Energy"),""));

  values->push_back
    (G4AttValue("Pos",G4BestUnit(finalPos,"Length"),""));

  return values;
}

void G4CMPElectrodeHit::Print() {
  G4cout << "  time " << finalTime/ns << " (nsec) : at " << finalPos
         << "  -- edep = " << EDep/eV << " [eV]" << G4endl;
}



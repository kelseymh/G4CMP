/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// 20200510  M. Kelsey -- G4CMP-201: Allocator must be thread-local
// 20220222  G4CMP-289 -- Thread-local allocator must be a pointer.

#ifndef G4CMPElectrodeHit_h
#define G4CMPElectrodeHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class G4AttDef;
class G4AttValue;

class G4CMPElectrodeHit : public G4VHit {
public:
  G4CMPElectrodeHit();
  int operator==(const G4CMPElectrodeHit &right) const;

  inline void *operator new(size_t);
  inline void operator delete(void *aHit);

  virtual void Draw();
  virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
  virtual std::vector<G4AttValue>* CreateAttValues() const;
  virtual void Print();

  void SetStartTime(G4double t) { startTime = t; }
  G4double GetStartTime() const { return startTime; }

  void SetFinalTime(G4double t) { finalTime = t; }
  G4double GetFinalTime() const { return finalTime; }

  void SetStartEnergy(G4double E) { startE = E; }
  G4double GetStartEnergy() const { return startE; }

  void SetEnergyDeposit(G4double E) { EDep = E; }
  G4double GetEnergyDeposit() const { return EDep; }

  void SetWeight(G4double w) { weight = w; }
  G4double GetWeight() const { return weight; }

  void SetStartPosition(G4ThreeVector xyz) { startPos = xyz; }
  G4ThreeVector GetStartPosition() const { return startPos; }

  void SetFinalPosition(G4ThreeVector xyz) { finalPos = xyz; }
  G4ThreeVector GetFinalPosition() const { return finalPos; }

  void SetTrackID(G4int id) { trackID = id; }
  G4int GetTrackID() const { return trackID; }

  void SetParticleName(G4String name) { particleName = name; }
  G4String GetParticleName() const { return particleName; }

private:
  G4double startTime;
  G4double finalTime;
  G4double startE;
  G4double EDep;
  G4double weight;
  G4ThreeVector startPos;
  G4ThreeVector finalPos;
  G4int trackID;
  G4String particleName;
};

typedef G4THitsCollection<G4CMPElectrodeHit> G4CMPElectrodeHitsCollection;

extern G4ThreadLocal G4Allocator<G4CMPElectrodeHit>* G4CMPElectrodeHitAllocator;

inline void* G4CMPElectrodeHit::operator new(size_t) {
  if (!G4CMPElectrodeHitAllocator)
    G4CMPElectrodeHitAllocator = new G4Allocator<G4CMPElectrodeHit>;
  return (void*)G4CMPElectrodeHitAllocator->MallocSingle();
}

inline void G4CMPElectrodeHit::operator delete(void* aHit) {
  G4CMPElectrodeHitAllocator->FreeSingle((G4CMPElectrodeHit*) aHit);
}

#endif	/* G4CMPElectrodeHit_h */

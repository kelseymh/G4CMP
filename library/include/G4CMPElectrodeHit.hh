/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

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
      ~G4CMPElectrodeHit();
      int operator==(const G4CMPElectrodeHit &right) const;

      inline void *operator new(size_t);
      inline void operator delete(void *aHit);

      virtual void Draw();
      virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
      virtual std::vector<G4AttValue>* CreateAttValues() const;
      virtual void Print();

  private:
    G4int trackID;
    G4String particleName;
    G4double startTime;
    G4double finalTime;
    G4double startE;
    G4double EDep;
    G4ThreeVector startPos;
    G4ThreeVector finalPos;

  public:
    inline void SetTrackID(G4int id) { trackID = id; }
    inline G4int GetTrackID() const { return trackID; }

    inline void SetParticleName(G4String name) { particleName = name; }
    inline G4String GetParticleName() const { return particleName; }

    inline void SetStartTime(G4double t) { startTime = t; }
    inline G4double GetStartTime() const { return startTime; }

    inline void SetFinalTime(G4double t) { finalTime = t; }
    inline G4double GetFinalTime() const { return finalTime; }

    inline void SetStartEnergy(G4double E) { startE = E; }
    inline G4double GetStartEnergy() const { return startE; }

    inline void SetEnergyDeposit(G4double E) { EDep = E; }
    inline G4double GetEnergyDeposit() const { return EDep; }

    inline void SetStartPosition(G4ThreeVector xyz) { startPos = xyz; }
    inline G4ThreeVector GetStartPosition() const { return startPos; }

    inline void SetFinalPosition(G4ThreeVector xyz) { finalPos = xyz; }
    inline G4ThreeVector GetFinalPosition() const { return finalPos; }
};

typedef G4THitsCollection<G4CMPElectrodeHit> G4CMPElectrodeHitsCollection;

extern G4Allocator<G4CMPElectrodeHit> G4CMPElectrodeHitAllocator;

inline void* G4CMPElectrodeHit::operator new(size_t) {
  void* aHit;
  aHit = (void*)G4CMPElectrodeHitAllocator.MallocSingle();
  return aHit;
}

inline void G4CMPElectrodeHit::operator delete(void* aHit) {
  G4CMPElectrodeHitAllocator.FreeSingle((G4CMPElectrodeHit*) aHit);
}

#endif	/* G4CMPElectrodeHit_h */

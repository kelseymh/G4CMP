#ifndef G4CMPElectrodeHit_h
#define G4CMPElectrodeHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class G4AttDef;
class G4AttValue;

class G4CMPElectrodeHit : public G4VHit
{
  public:

      G4CMPElectrodeHit();
      virtual ~G4CMPElectrodeHit();
      G4CMPElectrodeHit(const G4CMPElectrodeHit &right);
      const G4CMPElectrodeHit& operator=(const G4CMPElectrodeHit &right);
      int operator==(const G4CMPElectrodeHit &right) const;

      inline void *operator new(size_t);
      inline void operator delete(void *aHit);

      virtual void Draw();
      virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
      virtual std::vector<G4AttValue>* CreateAttValues() const;
      virtual void Print();

  private:
    G4int trackID;
    G4double charge;
    G4double trackTime;
    G4double startE;
    G4double EDep;
    G4ThreeVector startPos;
    G4ThreeVector finalPos;

  public:
    inline void SetTrackID(G4int id) { trackID = id; }
    inline G4int GetTrackID() const { return trackID; }

    inline void SetCharge(G4int q) { charge = q; }
    inline G4double GetCharge() const { return charge; }

    inline void SetStartEnergy(G4double E) { EDep = E; }
    inline G4double GetStartEnergy() const { return EDep; }

    inline void SetTrackTime(G4double t) { trackTime = t; }
    inline G4double GetTrackTime() const { return trackTime; }

    inline void SetEnergyDeposit(G4double E) { EDep = E; }
    inline G4double GetEnergyDeposit() const { return EDep; }

    inline void SetStartPosition(G4ThreeVector xyz) { startPos = xyz; }
    inline G4ThreeVector GetStartPosition() const { return startPos; }

    inline void SetFinalPosition(G4ThreeVector xyz) { finalPos = xyz; }
    inline G4ThreeVector GetFinalPosition() const { return finalPos; }
};

typedef G4THitsCollection<G4CMPElectrodeHit> G4CMPElectrodeHitsCollection;

extern G4Allocator<G4CMPElectrodeHit> G4CMPElectrodeHitAllocator;

inline void* G4CMPElectrodeHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*)G4CMPElectrodeHitAllocator.MallocSingle();
  return aHit;
}

inline void G4CMPElectrodeHit::operator delete(void* aHit)
{
  G4CMPElectrodeHitAllocator.FreeSingle((G4CMPElectrodeHit*) aHit);
}

#endif



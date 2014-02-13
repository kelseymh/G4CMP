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
      G4double time;
      G4double edep;
      G4double charge;    //Need charge for FETSim 9/27/12
      G4ThreeVector localPos;
      G4ThreeVector worldPos;

  public:
      inline void SetTime(G4double t) { time = t; }
      inline G4double GetTime() const { return time; }
      inline void SetEDep(G4double e) { edep = e; }
      inline G4double GetEDep() const { return edep; }
      inline void SetLocalPos(G4ThreeVector xyz) { localPos = xyz; }
      inline G4ThreeVector GetLocalPos() const { return localPos; }
      inline void SetWorldPos(G4ThreeVector xyz) { worldPos = xyz; }
      inline G4ThreeVector GetWorldPos() const { return worldPos; }
      inline void SetCharge(G4double c) { charge = c; } //Added by Rob for FETSim 9/27/12
      inline G4double GetCharge() { return charge; } //Added by Rob for FETSim 9/27/12
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



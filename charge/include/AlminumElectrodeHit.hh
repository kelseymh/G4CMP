#ifndef AlminumElectrodeHit_h
#define AlminumElectrodeHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class G4AttDef;
class G4AttValue;

class AlminumElectrodeHit : public G4VHit
{
  public:

      AlminumElectrodeHit();
      virtual ~AlminumElectrodeHit();
      AlminumElectrodeHit(const AlminumElectrodeHit &right);
      const AlminumElectrodeHit& operator=(const AlminumElectrodeHit &right);
      int operator==(const AlminumElectrodeHit &right) const;

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

typedef G4THitsCollection<AlminumElectrodeHit> AlminumElectrodeHitsCollection;

extern G4Allocator<AlminumElectrodeHit> AlminumElectrodeHitAllocator;

inline void* AlminumElectrodeHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*)AlminumElectrodeHitAllocator.MallocSingle();
  return aHit;
}

inline void AlminumElectrodeHit::operator delete(void* aHit)
{
  AlminumElectrodeHitAllocator.FreeSingle((AlminumElectrodeHit*) aHit);
}

#endif



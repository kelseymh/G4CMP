#ifndef PhononTrackInformation_h
#define PhononTrackInformation_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4Allocator.hh"
#include "G4VUserTrackInformation.hh"

class PhononTrackInformation : public G4VUserTrackInformation
{

public:
  PhononTrackInformation();
  PhononTrackInformation(const G4Track* aTrack);
  PhononTrackInformation(const PhononTrackInformation* aTrackInfo);
  PhononTrackInformation(G4ThreeVector kNew);
  virtual ~PhononTrackInformation();
  
  inline void *operator new(size_t);
  inline void operator delete(void *aTrackInfo);
  inline int operator == (const PhononTrackInformation& right) const {return (this==&right);}

  PhononTrackInformation& operator = (const PhononTrackInformation& right);

  //void SetSourceTrackInformation(const G4Track* aTrack);
  void Print() const;

private:
  G4ThreeVector k;

public:
  inline G4ThreeVector getK() const {return k;}
  inline void setK(G4ThreeVector kNew){k=G4ThreeVector(kNew);}

};

extern G4Allocator<PhononTrackInformation> pTrackInformationAllocator;

inline void* PhononTrackInformation::operator new(size_t)
{ void* aTrackInfo;
  aTrackInfo = (void*)pTrackInformationAllocator.MallocSingle();
  return aTrackInfo;

}

inline void PhononTrackInformation::operator delete(void *aTrackInfo)
{ pTrackInformationAllocator.FreeSingle((PhononTrackInformation*)aTrackInfo);}

#endif

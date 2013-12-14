#ifndef DriftingElectronTrackInformation_h
#define DriftingElectronTrackInformation_h 1

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4Allocator.hh"
#include "G4VUserTrackInformation.hh"

class DriftingElectronTrackInformation : public G4VUserTrackInformation{

public:
  DriftingElectronTrackInformation();
  DriftingElectronTrackInformation(int v);
  DriftingElectronTrackInformation(const DriftingElectronTrackInformation* aTrackInfo);
  virtual ~DriftingElectronTrackInformation();
  
  inline int operator == (const DriftingElectronTrackInformation& right) const {return(this==&right);}
  inline void* operator new(size_t);
  inline void operator delete(void *aTrackInfo);

  void Print() const;

private:
  int valley;

public:
  inline int getValley() const {return valley;}
  inline void setValley(int v){valley=v;}

};

extern G4Allocator<DriftingElectronTrackInformation> eTrackInformationAllocator;

inline void* DriftingElectronTrackInformation::operator new(size_t){
  void* aTrackInfo;
  aTrackInfo=(void*)eTrackInformationAllocator.MallocSingle();
  return aTrackInfo;
}

inline void DriftingElectronTrackInformation::operator delete(void *aTrackInfo){
  eTrackInformationAllocator.FreeSingle((DriftingElectronTrackInformation*)aTrackInfo);}

#endif

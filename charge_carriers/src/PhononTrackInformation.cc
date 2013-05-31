#include "PhononTrackInformation.hh"
#include "G4ios.hh"

G4Allocator<PhononTrackInformation> pTrackInformationAllocator;

PhononTrackInformation::PhononTrackInformation()
{
  k = G4ThreeVector(0., 0., 0.);
}

PhononTrackInformation::PhononTrackInformation(const G4Track* aTrack)
{
  k = G4ThreeVector(0.,0.,0.);
}

PhononTrackInformation::PhononTrackInformation(const PhononTrackInformation* aTrackInfo)
{
  k = aTrackInfo->getK(); 
}

PhononTrackInformation::PhononTrackInformation(G4ThreeVector kNew){
  k=kNew;
}

PhononTrackInformation::~PhononTrackInformation()
{ ; }

PhononTrackInformation& PhononTrackInformation::operator=(const PhononTrackInformation& aTrackInfo)
{
  k=aTrackInfo.k;

  return *this;

}

void PhononTrackInformation::Print() const
{
  G4cout<<"\nPhononTrackInformation::Print: Phonon k-vector is "<<k;
}

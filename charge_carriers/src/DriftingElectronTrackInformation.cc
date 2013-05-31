#include "DriftingElectronTrackInformation.hh"
#include "G4ios.hh"

G4Allocator<DriftingElectronTrackInformation> eTrackInformationAllocator;

DriftingElectronTrackInformation::DriftingElectronTrackInformation()
{
  valley = 1;
}

DriftingElectronTrackInformation::DriftingElectronTrackInformation(int v)
{
  valley = v;
}

DriftingElectronTrackInformation::DriftingElectronTrackInformation(const DriftingElectronTrackInformation* aTrackInfo)
{
  valley = aTrackInfo->getValley();
}

DriftingElectronTrackInformation::~DriftingElectronTrackInformation()
{ ; }

void DriftingElectronTrackInformation::Print() const{
  G4cout<<"\nDriftingElectronTrackInformation::Print: Valley: "<<valley;
}

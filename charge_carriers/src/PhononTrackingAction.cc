
#include "PhononTrackingAction.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "PhononTrackInformation.hh"

void PhononTrackingAction::PreUserTrackingAction(const G4Track* aTrack){
  ;
}

void PhononTrackingAction::PostUserTrackingAction(const G4Track* aTrack){

  G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();

  if(secondaries)
    {
      PhononTrackInformation* info = (PhononTrackInformation*)(aTrack->GetUserInformation());
      //      size_t nSeco = secondaries->size();
      //      if(nSeco>0){
	//	for(size_t i=0;i<nSeco;i++){
	  //	  PhononTrackInformation* infoNew = new PhononTrackInformation(info);
	  //	  //	  G4cout<<"\nSetting secondary number "<<i<<" info pointer address"<<infoNew;
	  //(*secondaries)[i]->SetUserInformation(infoNew);
	  //}
	//      }
    }


}

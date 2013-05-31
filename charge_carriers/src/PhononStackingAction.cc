
#include "PhononStackingAction.hh"
#include "PhononTrackInformation.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "LatticeManager2.hh"
#include "G4RandomDirection.hh"
#include "G4ThreeVector.hh"
#include "TPhononFast.hh"
#include "TPhononSlow.hh"
#include "LPhonon.hh"

PhononStackingAction::PhononStackingAction()
{;}

PhononStackingAction::~PhononStackingAction()
{;}

G4ClassificationOfNewTrack PhononStackingAction::ClassifyNewTrack(const G4Track* aTrack)
{
  G4ClassificationOfNewTrack classification = fUrgent;

  //PhononTrackInformation* info;
  
  if(aTrack->GetParentID() == 0)
    {
      //info = new PhononTrackInformation();
      G4Track* theTrack=(G4Track*)aTrack;
      G4ThreeVector Ran = G4RandomDirection();
      G4int pol=0;
      if(theTrack->GetDefinition()==LPhonon::PhononDefinition()) pol = 0;
      else if(theTrack->GetDefinition()==TPhononSlow::PhononDefinition()) pol = 1;
      else if(theTrack->GetDefinition()==TPhononFast::PhononDefinition()) pol = 2;

      theTrack->SetUserInformation(new PhononTrackInformation(Ran));     
      theTrack->SetMomentumDirection(LatticeManager2::mapKtoVDir(theTrack->GetVolume(),pol,Ran));
      theTrack->SetMomentumDirection(G4RandomDirection());
    }
  return classification;
}

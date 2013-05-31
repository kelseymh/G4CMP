#include "XtalFieldManager.hh"
#include "DriftingElectronTrackInformation.hh"

#include "G4ParticleDefinition.hh"
#include "G4Track.hh"



    


void XtalFieldManager::ConfigureForTrack(const G4Track* aTrack){

  if(aTrack->GetDefinition()->GetParticleName()=="DriftingElectron"){
    int valley = ((DriftingElectronTrackInformation*) aTrack->GetUserInformation())->getValley();
    switch(valley){
      case 1:
	SetChordFinder(valley1ChordFinder);
	break;

    case 2:
      SetChordFinder(valley2ChordFinder);
      break;

    case 3:
      SetChordFinder(valley3ChordFinder);
      break;
 
    case 4:
      SetChordFinder(valley4ChordFinder);
      break;
    }

  } else {
    SetChordFinder(normalChordFinder);
  }
  

}



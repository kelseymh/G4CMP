/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#include "ChannelingTrackingAction.hh"

#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4TrackVector.hh"

#include "ChannelingParticleUserInfo.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChannelingTrackingAction::PreUserChannelingTrackingAction(const G4Track* aTrack){
    ChannelingParticleUserInfo *  trackInfo( static_cast< ChannelingParticleUserInfo * >(aTrack->GetUserInformation() ) );
    
    if ( trackInfo ){
        return;
    }
    else{
        G4Track *  theTrack( const_cast< G4Track * >( aTrack ) );
        trackInfo = new ChannelingParticleUserInfo();
        theTrack->SetUserInformation( trackInfo );
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ChannelingTrackingAction::PostUserChannelingTrackingAction(const G4Track*){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

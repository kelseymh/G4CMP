
#include "DriftingElectronStackingAction.hh"
#include "DriftingElectronTrackInformation.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "DriftingElectron.hh"
#include "Randomize.hh"

DriftingElectronStackingAction::DriftingElectronStackingAction()
{ ; }

DriftingElectronStackingAction::~DriftingElectronStackingAction()
{ ; }

G4ClassificationOfNewTrack DriftingElectronStackingAction::ClassifyNewTrack(const G4Track* aTrack)
{

  G4ClassificationOfNewTrack classification = fUrgent;

  G4Track* theTrack = (G4Track*)aTrack;

  int valley = (int) (G4UniformRand()*4 + 1.0);
  
  theTrack->SetUserInformation(new DriftingElectronTrackInformation(valley));

  return classification;

}

#ifndef DriftingElectronStackingAction_h
#define DriftingElectronStackingAction_h 1

#include "globals.hh"
#include "G4UserStackingAction.hh"

class G4Track;

class DriftingElectronStackingAction : public G4UserStackingAction
{

public:
  DriftingElectronStackingAction();
  virtual ~DriftingElectronStackingAction();

public:
  virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);

};


#endif

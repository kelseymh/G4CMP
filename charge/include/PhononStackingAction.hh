#ifndef PhononStackingAction_h
#define PhononStackingAction_h 1

#include "globals.hh"
#include "G4UserStackingAction.hh"

class G4Track;

class PhononStackingAction : public G4UserStackingAction
{
public:
  PhononStackingAction();
  virtual ~PhononStackingAction();

public:
  virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);


};





#endif

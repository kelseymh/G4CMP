#ifndef PhononTrackingAction_h
#define PhononTrackingAction_h 1

#include "G4UserTrackingAction.hh"

class PhononTrackingAction : public G4UserTrackingAction{

public:
  PhononTrackingAction(){};
  virtual ~PhononTrackingAction(){};
  
  virtual void PreUserTrackingAction(const G4Track*);
  virtual void PostUserTrackingAction(const G4Track*);

};

#endif

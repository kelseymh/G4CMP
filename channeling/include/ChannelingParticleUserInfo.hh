#ifndef ChannelingUserInfo_h
#define ChannelingUserInfo_h 1

#include "globals.hh"
#include "G4VUserTrackInformation.hh"

class ChannelingParticleUserInfo : public G4VUserTrackInformation
{

public:

  ChannelingParticleUserInfo();
  ~ChannelingParticleUserInfo();

  void SetChanneling(bool flag);
  bool GetChanneling();

  void SetChannelingFactor(G4double);
  G4double GetChannelingFactor();
  G4double GetPreStepChannelingFactor();

private:

  G4double channelingFactor;
  G4double preStepChannelingFactor;
  bool channelingFlag;

};



#endif

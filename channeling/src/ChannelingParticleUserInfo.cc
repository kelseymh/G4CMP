#include "ChannelingParticleUserInfo.hh"

ChannelingParticleUserInfo::ChannelingParticleUserInfo()
{
  channelingFlag = false;
  channelingFactor = 1.0;
  preStepChannelingFactor = 1.0;
}

ChannelingParticleUserInfo::~ChannelingParticleUserInfo()
{
  channelingFlag = false;
  channelingFactor = 1.0;
}

void ChannelingParticleUserInfo::SetChanneling(bool flag)
{
  channelingFlag = flag;
}

void ChannelingParticleUserInfo::SetChannelingFactor(G4double factor)
{
  preStepChannelingFactor = channelingFactor;
  channelingFactor = factor;
}

G4double ChannelingParticleUserInfo::GetChannelingFactor()
{
  return channelingFactor;
}

G4double ChannelingParticleUserInfo::GetPreStepChannelingFactor()
{
  return preStepChannelingFactor;
}

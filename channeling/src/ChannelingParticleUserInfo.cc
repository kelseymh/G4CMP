#include "ChannelingParticleUserInfo.hh"

ChannelingParticleUserInfo::ChannelingParticleUserInfo()
{
  channelingFlag = false;
  channelingFactor = 1.0;
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
  channelingFactor = factor;
}

G4double ChannelingParticleUserInfo::GetChannelingFactor()
{
  return channelingFactor;
}

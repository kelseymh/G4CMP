/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file electromagnetic/TestEm5/src/ChannelingStepLimiter.cc
/// \brief Implementation of the ChannelingStepLimiter class
//
// $Id: 6ed964dec2b4209692f9dc41bc077ea0d2cd65fe $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ChannelingStepLimiter.hh"
#include "ChannelingStepMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ChannelingStepLimiter::ChannelingStepLimiter(const G4String& processName)
 : G4VDiscreteProcess(processName),fMaxChargedStep(DBL_MAX)
{
  fMess = new ChannelingStepMessenger(this);
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ChannelingStepLimiter::~ChannelingStepLimiter() { delete fMess; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool ChannelingStepLimiter::IsApplicable(const G4ParticleDefinition& particle) 
{ 
  return (particle.GetPDGCharge() != 0.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    
void ChannelingStepLimiter::SetMaxStep(G4double step) { fMaxChargedStep = step;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double ChannelingStepLimiter::PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
                                                  G4double,
                                                  G4ForceCondition* condition )
{
  // condition is set to "Not Forced"
  *condition = NotForced;
  
  G4double ProposedStep = DBL_MAX;

  if((fMaxChargedStep > 0.) &&
     (aTrack.GetVolume() != 0) &&
     (aTrack.GetVolume()->GetName() == "Target"))
     ProposedStep = fMaxChargedStep;

  return ProposedStep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VParticleChange* ChannelingStepLimiter::PostStepDoIt(const G4Track& aTrack, const G4Step&)
{
   // do nothing
   aParticleChange.Initialize(aTrack);
   return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



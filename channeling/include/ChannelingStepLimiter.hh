/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file electromagnetic/TestEm5/include/ChannelingStepLimiter.hh
/// \brief Definition of the ChannelingStepLimiter class
//
// $Id: b0d7672042db3952ca9b924d17454d3f5db3a115 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef ChannelingStepLimiter_h
#define ChannelingStepLimiter_h 1

#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "G4ParticleDefinition.hh"
#include "G4Step.hh"

class ChannelingStepMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ChannelingStepLimiter : public G4VDiscreteProcess
{
  public:     

     ChannelingStepLimiter(const G4String& processName ="UserChannelingStepLimiter");
    ~ChannelingStepLimiter();

     virtual G4bool   IsApplicable(const G4ParticleDefinition&);    
     void     SetMaxStep(G4double);
     G4double GetMaxStep() {return fMaxChargedStep;};
     
     virtual G4double PostStepGetPhysicalInteractionLength(const G4Track& track,
                                             G4double   previousStepSize,
                                             G4ForceCondition* condition);

     virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

     virtual G4double GetMeanFreePath(const G4Track&,G4double,G4ForceCondition*)
       {return 0.;};     // it is not needed here !

  private:

     G4double    fMaxChargedStep;
     ChannelingStepMessenger* fMess;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


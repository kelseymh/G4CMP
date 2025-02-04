/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPBogoliubovQPDiffusionTimeeStepperProcess.cc
/// \brief Implementation of the G4CMPQPDiffusionTimeStepperProcess class
//
// $Id$
//


#include "G4CMPQPDiffusionTimeStepperProcess.hh"
#include "G4CMPQPDiffusionTimeStepperRate.hh"
#include "G4CMPSCUtils.hh"
#include "G4BogoliubovQP.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VParticleChange.hh"
#include "G4RandomDirection.hh"
#include "G4CMPUtils.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4CMPSecondaryUtils.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Constructor and destructor 
G4CMPQPDiffusionTimeStepperProcess::G4CMPQPDiffusionTimeStepperProcess(const G4String& aName)
  : G4VBogoliubovQPProcess(aName,fQPDiffusionTimeStepper)
{
  G4cout << "REL -- In QPDiffusionTimeStepper::Constructor()" << G4endl;
  UseRateModel(new G4CMPQPDiffusionTimeStepperRate);
}

G4CMPQPDiffusionTimeStepperProcess::~G4CMPQPDiffusionTimeStepperProcess()
{
  G4cout << "REL -- In QPDiffusionTimeStepper::Destructor()" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4CMPQPDiffusionTimeStepperProcess::SetVerboseLevel(G4int vb) {
  verboseLevel = vb;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4VParticleChange* G4CMPQPDiffusionTimeStepperProcess::PostStepDoIt(const G4Track& aTrack,
								       const G4Step& aStep)
{
  aParticleChange.Initialize(aTrack);

  //Pseudocode
  //1. Determine if we're on a boundary surface. If we are, kill the event -- if the code is working properly, this should
  //   basically never happen... REL May need to revisit this. I think it's still true but will have to check against the QP transport stuff.
  G4StepPoint * postStepPoint = aStep.GetPostStepPoint();
  if (postStepPoint->GetStepStatus() == fGeomBoundary ||
      postStepPoint->GetStepStatus() == fWorldBoundary) {
    G4ExceptionDescription msg;
    msg << "For some reason we're running post-step do it for the QPDiffusionTimeStepper process and we find ourselves on a boundary. Should this ever happen?";
    G4Exception("G4CMPQPDiffusionTimeStepperProcess::PostStepDoIt", "QPDiffusionTimeStepper001",EventMustBeAborted,msg);
    return &aParticleChange;		
  }
   
  //2. Identify the current QP's energy and velocity. No changes made here.
  G4double qpEnergy = GetKineticEnergy(aTrack);
  G4double velocity = aStep.GetPostStepPoint()->GetVelocity();
  aParticleChange.ProposeEnergy(qpEnergy);
  aParticleChange.ProposeVelocity(velocity);
  
  //3. Randomize the final state momentum of the QP in XY.
  RandomizeFinalStateMomentumDirectionInXY();
  
  //4. Do the clear interaction lengths thing because we do still have a particle here.
  ClearNumberOfInteractionLengthLeft();		// All processes should do this!
  
  //5. Return the particle change
  return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4bool G4CMPQPDiffusionTimeStepperProcess::IsApplicable(const G4ParticleDefinition& aPD) {
  G4cout << "REL -- In QPDiffusionTimeStepperProcess::IsApplicable()" << G4endl;
  // Allow all phonon types, because type is changed during tracking
  return G4VBogoliubovQPProcess::IsApplicable(aPD);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//Pass-through to G4CMPVProcess class
G4double G4CMPQPDiffusionTimeStepperProcess::GetMeanFreePath(const G4Track& trk, G4double prevstep, G4ForceCondition* cond)
{
  return G4CMPVProcess::GetMeanFreePath(trk,prevstep,cond);
}

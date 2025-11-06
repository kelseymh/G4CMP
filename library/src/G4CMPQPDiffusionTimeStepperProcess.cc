/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPQPDiffusionTimeeStepperProcess.cc
/// \brief Implementation of the G4CMPQPDiffusionTimeStepperProcess class
//
// $Id$
//

#include "G4CMPQPDiffusionTimeStepperProcess.hh"
#include "G4CMPBogoliubovQP.hh"
#include "G4CMPQPDiffusionTimeStepperRate.hh"
#include "G4CMPSCUtils.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VParticleChange.hh"
#include "G4RandomDirection.hh"
#include "G4CMPUtils.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4CMPSecondaryUtils.hh"

// Constructor and destructor 
G4CMPQPDiffusionTimeStepperProcess::
G4CMPQPDiffusionTimeStepperProcess(const G4String& aName)
  : G4VQPProcess(aName,fQPDiffusionTimeStepper) {
  UseRateModel(new G4CMPQPDiffusionTimeStepperRate);
}

G4CMPQPDiffusionTimeStepperProcess::~G4CMPQPDiffusionTimeStepperProcess() {
}

void G4CMPQPDiffusionTimeStepperProcess::SetVerboseLevel(G4int vb) {
  verboseLevel = vb;
}

G4VParticleChange* G4CMPQPDiffusionTimeStepperProcess::
PostStepDoIt(const G4Track& aTrack,const G4Step& aStep) {

  aParticleChange.Initialize(aTrack);

  //Pseudocode
  //1. Determine if we're on a boundary surface. If we are, kill the event --
  //   if the code is working properly, this should basically never happen...
  //   May need to revisit this. I think it's still true but will have to check
  //   against the QP transport stuff.
  G4StepPoint * postStepPoint = aStep.GetPostStepPoint();
  if (postStepPoint->GetStepStatus() == fGeomBoundary ||
      postStepPoint->GetStepStatus() == fWorldBoundary) {
    G4ExceptionDescription msg;
    msg << "For some reason we're running post-step do it for the "
	<< "QPDiffusionTimeStepper process and we find ourselves on a boundary."
	<< "Should this ever happen?";
    G4Exception("G4CMPQPDiffusionTimeStepperProcess::PostStepDoIt",
		"QPDiffusionTimeStepper001",EventMustBeAborted,msg);
    return &aParticleChange;		
  }
   
  //2. Identify the current QP's energy and velocity. No changes made here.
  G4double qpEnergy = GetKineticEnergy(aTrack);
  G4double velocity = aStep.GetPostStepPoint()->GetVelocity();
  aParticleChange.ProposeEnergy(qpEnergy);
  aParticleChange.ProposeVelocity(velocity);
  
  //3. Randomize the final state momentum of the QP in XY.
  RandomizeFinalStateMomentumDirectionInXY();
  
  //4. Do the clear interaction lengths thing because we do still have a
  //particle here.
  ClearNumberOfInteractionLengthLeft();	// All processes should do this!
  
  //5. Return the particle change
  return &aParticleChange;
}

G4bool G4CMPQPDiffusionTimeStepperProcess::
IsApplicable(const G4ParticleDefinition& aPD) {   
  return G4VQPProcess::IsApplicable(aPD);
}

//Pass-through to G4CMPVProcess class
G4double G4CMPQPDiffusionTimeStepperProcess::
GetMeanFreePath(const G4Track& trk,G4double prevstep,G4ForceCondition* cond) {
  return G4CMPVProcess::GetMeanFreePath(trk,prevstep,cond);
}

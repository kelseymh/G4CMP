/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPQPLocalTrappingProcess.cc
/// \brief Implementation of the G4CMPQPLocalTrappingProcess class
//
//

#include "G4CMPQPLocalTrappingProcess.hh"
#include "G4CMPQPLocalTrappingRate.hh"
#include "G4CMPSCUtils.hh"
#include "G4QP.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VParticleChange.hh"
#include "G4RandomDirection.hh"
#include "G4CMPUtils.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4CMPSecondaryUtils.hh"



// Constructor and destructor 
G4CMPQPLocalTrappingProcess::G4CMPQPLocalTrappingProcess(const G4String& aName)
  : G4VQPProcess(aName,fQPLocalTrappingProcess) {  
  UseRateModel(new G4CMPQPLocalTrappingRate);
}

G4CMPQPLocalTrappingProcess::~G4CMPQPLocalTrappingProcess() {
}


G4VParticleChange* G4CMPQPLocalTrappingProcess::
PostStepDoIt(const G4Track& aTrack,const G4Step& aStep) {

  //Debugging
  if (verboseLevel > 5) {
    G4cout << "-- G4CMPQPLocalTrappingProcess::PostStepDoIt --" << G4endl;
  }
  
  aParticleChange.Initialize(aTrack);
  
  G4StepPoint* postStepPoint = aStep.GetPostStepPoint();
  if (postStepPoint->GetStepStatus()==fGeomBoundary ||
      postStepPoint->GetStepStatus()==fWorldBoundary) {
    return &aParticleChange; // Don't want to reset IL
  }

  if (verboseLevel) G4cout << GetProcessName() << "::PostStepDoIt" << G4endl;
  if (verboseLevel>1) {
    G4StepPoint* preStepPoint = aStep.GetPreStepPoint();
    G4cout << " Track " << aTrack.GetDefinition()->GetParticleName()
	   << " vol " << aTrack.GetTouchable()->GetVolume()->GetName()
	   << " prePV " << preStepPoint->GetPhysicalVolume()->GetName()
	   << " postPV " << postStepPoint->GetPhysicalVolume()->GetName()
	   << " step-length " << aStep.GetStepLength()
	   << G4endl;
  }

  //Since this local trapping, this is very simple -- we just kill the track.
  aParticleChange.ProposeTrackStatus(fStopAndKill);

  //Do the clear interaction lengths thing because we do still have a particle
  //here.
  ClearNumberOfInteractionLengthLeft();

  //5. Return the particle change
  return &aParticleChange;
}

//Pass-through to G4CMPVProcess class
G4double G4CMPQPLocalTrappingProcess::
GetMeanFreePath(const G4Track& trk,G4double prevstep,G4ForceCondition* cond) {
  
  return G4CMPVProcess::GetMeanFreePath(trk,prevstep,cond);
}

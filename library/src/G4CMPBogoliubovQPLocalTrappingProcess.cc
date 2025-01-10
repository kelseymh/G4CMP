/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPBogoliubovQPLocalTrappingProcess.cc
/// \brief Implementation of the G4CMPBogoliubovQPLocalTrappingProcess class
//
//

#include "G4CMPBogoliubovQPLocalTrappingProcess.hh"
#include "G4CMPBogoliubovQPLocalTrappingRate.hh"
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
G4CMPBogoliubovQPLocalTrappingProcess::G4CMPBogoliubovQPLocalTrappingProcess(const G4String& aName)
  : G4VBogoliubovQPProcess(aName,fQPLocalTrappingProcess)
{  
  UseRateModel(new G4CMPBogoliubovQPLocalTrappingRate);
}

G4CMPBogoliubovQPLocalTrappingProcess::~G4CMPBogoliubovQPLocalTrappingProcess()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4VParticleChange* G4CMPBogoliubovQPLocalTrappingProcess::PostStepDoIt(const G4Track& aTrack,
								       const G4Step& aStep)
{
  G4cout << "REL in BogoliubovQPLocalTrapping process poststepdoit." << G4endl;
  
  aParticleChange.Initialize(aTrack);

  
  G4StepPoint* postStepPoint = aStep.GetPostStepPoint();
  if (postStepPoint->GetStepStatus()==fGeomBoundary ||
      postStepPoint->GetStepStatus()==fWorldBoundary) {
    return &aParticleChange;			// Don't want to reset IL
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

  //Do the clear interaction lengths thing because we do still have a particle here.
  ClearNumberOfInteractionLengthLeft();

  //5. Return the particle change
  return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//Pass-through to G4CMPVProcess class
G4double G4CMPBogoliubovQPLocalTrappingProcess::GetMeanFreePath(const G4Track& trk, G4double prevstep, G4ForceCondition* cond)
{
  //Need this to come first, so that it actually attempts a superconductor update.
  G4double mfpBase = G4CMPVProcess::GetMeanFreePath(trk,prevstep,cond);

  //Here, before we try to run this, check to see if all of the relevant crystal parameters are defined. If they aren't,
  //throw an exception.
  G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();
  G4LatticePhysical* theLat;
  G4VPhysicalVolume* volume = trk.GetVolume();
  theLat = LM->GetLattice(volume);
  G4double localTrappingMFP = theLat->GetSCQPLocalTrappingMFP();
  if( localTrappingMFP == DBL_MAX ){
    G4ExceptionDescription msg;
    msg << "Noticed that in the mean free path calculation step for the BogoliubovQP local trapping, you have incorrectly defined or omitted the trapping MFP. In other words, you don't have enough input information in your config.txt file to run the trapping physics correctly. Please add the sc_qptrapmfp parameter with a physical value and unit to config.txt";
    G4Exception("G4CMPBogoliubovQPLocalTrapping::GetMeanFreePath", "BogoliubovQPLocalTrapping001",FatalException, msg);
  }

  //If we don't trigger that exception, continue.
  return mfpBase;
}

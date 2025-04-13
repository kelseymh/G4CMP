/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPParticleChangeForPhonon.cc
/// \brief Implementation of the G4CMPParticleChangeForPhonon class
//
// $Id$
//
// 20250410 Implement ParticleChange for phonons to handle displaced reflections

#include "G4CMPParticleChangeForPhonon.hh"
#include "G4TouchableHandle.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"

G4CMPParticleChangeForPhonon::G4CMPParticleChangeForPhonon() {;}

G4Step* G4CMPParticleChangeForPhonon::UpdateStepForPostStep(G4Step* pStep) {
  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint();
  // Update next volume if touchable has been proposed
  if (updateVol) {
    // Update next touchable
    pPostStepPoint->SetTouchableHandle(theTouchableHandle);
    pPostStepPoint->SetMaterial(theMaterialChange);
    pPostStepPoint->SetMaterialCutsCouple(theMaterialCutsCoupleChange);
    pPostStepPoint->SetSensitiveDetector(theSensitiveDetectorChange);

    if(this->GetFirstStepInVolume()) {
      pStep->SetFirstStepFlag();
    } else {
      pStep->ClearFirstStepFlag();
    }

    if(this->GetLastStepInVolume()) {
      pStep->SetLastStepFlag();
    } else {
      pStep->ClearLastStepFlag();
    }
    
    // Reset updateVol
    updateVol = false;
  }

  // Call base class function
  return G4ParticleChange::UpdateStepForPostStep(pStep);
}

void G4CMPParticleChangeForPhonon::ProposeTouchableHandle(G4TouchableHandle nextTouchableHandle) {
  // Set the new touchable handle and corresponding values
  theTouchableHandle = nextTouchableHandle;
  theMaterialChange = theTouchableHandle->GetVolume()->GetLogicalVolume()->GetMaterial();
  theSensitiveDetectorChange = theTouchableHandle->GetVolume()->GetLogicalVolume()->GetSensitiveDetector();
  updateVol = true;
}
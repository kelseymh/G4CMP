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
// 20250413 Add Initialize() implementation to reset updateVol flag

#include "G4CMPParticleChangeForPhonon.hh"
#include "G4TouchableHandle.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4VPhysicalVolume.hh"


// Copy operations should call back to base

G4CMPParticleChangeForPhonon::
G4CMPParticleChangeForPhonon(const G4CMPParticleChangeForPhonon& right)
  : G4ParticleChange(right), theTouchableHandle(right.theTouchableHandle),
    updateVol(right.updateVol) {;}

G4CMPParticleChangeForPhonon&
G4CMPParticleChangeForPhonon::operator=(const G4CMPParticleChangeForPhonon& right) {
  G4ParticleChange::operator=(right);
  updateVol = right.updateVol;
  theTouchableHandle = right.theTouchableHandle;
  return *this;
}


// Ensure that local flags are initialized at each step

void G4CMPParticleChangeForPhonon::Initialize(const G4Track& track) {
  updateVol = false;
  theTouchableHandle = track.GetTouchableHandle();
  
  G4ParticleChange::Initialize(track);
}


// Populate G4Step with results of process' PostStepDoIt()
  
G4Step* G4CMPParticleChangeForPhonon::UpdateStepForPostStep(G4Step* pStep) {
  if (updateVol) {    // Update next volume if touchable has been proposed
    G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint();
    G4LogicalVolume* LV = theTouchableHandle->GetVolume()->GetLogicalVolume();
    pPostStepPoint->SetTouchableHandle(theTouchableHandle);
    pPostStepPoint->SetMaterial(LV->GetMaterial());
    pPostStepPoint->SetMaterialCutsCouple(LV->GetMaterialCutsCouple());
    pPostStepPoint->SetSensitiveDetector(LV->GetSensitiveDetector());

    // Reset updateVol
    updateVol = false;
  }

  // Call base class function
  return G4ParticleChange::UpdateStepForPostStep(pStep);
}


// Include local information in printout

void G4CMPParticleChangeForPhonon::DumpInfo() const {
  G4ParticleChange::DumpInfo();

  G4cout << "        updateVol : " << updateVol << G4endl;
  if (updateVol) {
    G4cout << "        theTouchableHandle for PV "
	   << theTouchableHandle->GetVolume()->GetName() << G4endl;
  }
}

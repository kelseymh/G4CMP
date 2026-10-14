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
// 20251116 For G4 11, explicitly remove the copy operators to match base.
// 20251128 Discard touchable contents after updating.
// 20260212 G4CMP-585 Only discard touchable if it was modified (fix mem leak)
// 20260929 G4CMP-665: Remove G4TouchableHandle to avoid memory leak and seg
//          faults and rename function to UpdateTouchableHandle.

#include "G4CMPParticleChangeForPhonon.hh"
#include "G4VTouchable.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4VPhysicalVolume.hh"

// Ensure that local flags are initialized at each step

void G4CMPParticleChangeForPhonon::Initialize(const G4Track& track) {
  // Use base class's method first.
  G4ParticleChange::Initialize(track);
  updateVol = false;
}


// Populate G4Step with results of process' PostStepDoIt()
  
G4Step* G4CMPParticleChangeForPhonon::UpdateStepForPostStep(G4Step* pStep) {
  if (updateVol) {    // Update next volume if touchable has been proposed
    G4LogicalVolume* LV = pStep->GetPreStepPoint()->GetTouchableHandle()
                               ->GetVolume()->GetLogicalVolume();

    G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint();
    pPostStepPoint->SetTouchableHandle(pStep->GetPreStepPoint()->GetTouchableHandle());
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
}

/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPParticleChangeForPhonon.hh
/// \brief Definition of the G4CMPParticleChangeForPhonon class
///   Concrete class for ParticleChange for handling phonon reflections
//
// $Id$
//
// 20250410 Implement ParticleChange for phonons to handle displaced reflections

#ifndef G4CMPParticleChangeForPhonon_hh
#define G4CMPParticleChangeForPhonon_hh 1

#include "G4TouchableHandle.hh"
#include "G4ParticleChange.hh"

class G4CMPParticleChangeForPhonon final : public G4ParticleChange {
public:
  G4CMPParticleChangeForPhonon();
  ~G4CMPParticleChangeForPhonon() override = default;

  G4CMPParticleChangeForPhonon(const G4CMPParticleChangeForPhonon& right) = delete;
  G4CMPParticleChangeForPhonon& operator=(const G4CMPParticleChangeForPhonon& right) = delete;

  // --- Methods for updating G4Step ---
  G4Step* UpdateStepForPostStep(G4Step* pStep) final;
  
  // --- Methods for proposing PostStep volume ---
  void ProposeTouchableHandle(G4TouchableHandle nextTouchableHandle);

private:
  G4TouchableHandle theTouchableHandle;
  G4Material* theMaterialChange = nullptr;
  const G4MaterialCutsCouple* theMaterialCutsCoupleChange = nullptr;
  G4VSensitiveDetector* theSensitiveDetectorChange = nullptr;
  G4bool updateVol = false;
};

#endif /* G4CMPParticleChangeForPhonon_hh */
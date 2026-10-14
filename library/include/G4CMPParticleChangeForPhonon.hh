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
// 20250413 Add Initialize() implementation to reset updateVol flag, add
//          missing copy operations, may be needed
// 20251116 For G4 11, explicitly remove the copy operators to match base.
// 20251128 Implement empty destructor to avoid deleting G4TouchableHandle.
// 20260929 G4CMP-665: Remove G4TouchableHandle to avoid memory leak and seg
//          faults and rename function to UpdateTouchableHandle.

#ifndef G4CMPParticleChangeForPhonon_hh
#define G4CMPParticleChangeForPhonon_hh 1

#include "G4TouchableHandle.hh"
#include "G4ParticleChange.hh"


class G4CMPParticleChangeForPhonon final : public G4ParticleChange {
public:
  G4CMPParticleChangeForPhonon() : G4ParticleChange() {;}
  virtual ~G4CMPParticleChangeForPhonon() override = default;

  // ParticleChange cannot be copied
  G4CMPParticleChangeForPhonon(const G4CMPParticleChangeForPhonon& right) = delete;
  G4CMPParticleChangeForPhonon& operator=(const G4CMPParticleChangeForPhonon& right) = delete;

  // Ensure that local flags are cleared between steps
  virtual void Initialize(const G4Track& track) override;
  
  // --- Methods for updating G4Step ---
  G4Step* UpdateStepForPostStep(G4Step* pStep) final;
  
  // --- Methods for proposing PostStep volume ---
  void UpdateTouchableHandle() {
    updateVol = true;
  }

  // Include local information in printout
  virtual void DumpInfo() const override;
  
private:
  G4bool updateVol = false;		// Only set if touchable is changed
};

#endif /* G4CMPParticleChangeForPhonon_hh */

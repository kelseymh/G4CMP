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
//		missing copy operations, may be needed

#ifndef G4CMPParticleChangeForPhonon_hh
#define G4CMPParticleChangeForPhonon_hh 1

#include "G4TouchableHandle.hh"
#include "G4ParticleChange.hh"


class G4CMPParticleChangeForPhonon final : public G4ParticleChange {
public:
  G4CMPParticleChangeForPhonon() : G4ParticleChange() {;}
  virtual ~G4CMPParticleChangeForPhonon() override = default;

  // Ensure that local flags are cleared between steps
  virtual void Initialize(const G4Track& track) override;
  
  // --- Methods for updating G4Step ---
  G4Step* UpdateStepForPostStep(G4Step* pStep) final;
  
  // --- Methods for proposing PostStep volume ---
  void ProposeTouchableHandle(const G4TouchableHandle& nextTouchableHandle) {
    theTouchableHandle = nextTouchableHandle;
    updateVol = true;
  }
  
  const G4TouchableHandle& GetTouchableHandle() const {
    return theTouchableHandle;
  }

  // Include local information in printout
  virtual void DumpInfo() const override;
  
protected:	// Subclasses permitted to copy themselves
  G4CMPParticleChangeForPhonon(const G4CMPParticleChangeForPhonon& right);
  G4CMPParticleChangeForPhonon& operator=(const G4CMPParticleChangeForPhonon& right);

private:
  G4TouchableHandle theTouchableHandle;
  G4bool updateVol = false;		// Only set if touchable is changed
};

#endif /* G4CMPParticleChangeForPhonon_hh */

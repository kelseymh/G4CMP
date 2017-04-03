/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPUtils.hh
/// \brief Namespace for general purpose or static utilities supporting
///	G4CMP.  Functions not dependent on current track/step info will
///     be moved here from G4CMPProcessUtils.
//
// $Id$

#ifndef G4CMPUtils_hh
#define G4CMPUtils_hh 1

#include "G4ThreeVector.hh"
#include "globals.hh"

class G4CMPElectrodeHit;
class G4LatticePhysical;
class G4ParticleDefinition;
class G4Step;
class G4Track;


namespace G4CMP {
  template <class T> G4int sign(T val) {
    return (T(0) < val) - (val < T(0));
  }

  // Select phonon mode randomly from density of states
  G4int ChoosePhononPolarization(const G4LatticePhysical* lattice);
  G4int ChoosePhononPolarization(G4double Ldos, G4double STdos, G4double FTdos);

  // Randomly choose valley for charge carrier
  G4int ChooseValley(const G4LatticePhysical* lattice);

  // Identify G4CMP particle categories
  G4bool IsPhonon(const G4Track* track);
  G4bool IsElectron(const G4Track* track);
  G4bool IsHole(const G4Track* track);
  G4bool IsChargeCarrier(const G4Track* track);

  G4bool IsPhonon(const G4ParticleDefinition* pd);
  G4bool IsElectron(const G4ParticleDefinition* pd);
  G4bool IsHole(const G4ParticleDefinition* pd);
  G4bool IsChargeCarrier(const G4ParticleDefinition* pd);

  // Throw biasing decision for particle production and return weight
  G4double ChooseWeight(const G4ParticleDefinition* pd);
  G4double ChoosePhononWeight();
  G4double ChooseChargeWeight();

  // Create a Hit from a G4Step. Less error prone to use this helper.
  void FillHit(const G4Step*, G4CMPElectrodeHit*);

  // Phonons reflect difusively from surfaces.
  G4ThreeVector LambertReflection(const G4ThreeVector& surfNorm);

  // Test that a phonon's wave vector relates to an inward velocity.
  G4bool PhononVelocityIsInward(const G4LatticePhysical* lattice,
                                G4int polarization,
                                G4ThreeVector waveVector,
                                G4ThreeVector surfNorm);

}

#endif	/* G4CMPUtils_hh */

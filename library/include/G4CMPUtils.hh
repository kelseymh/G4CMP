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

#include "G4ThreeVector.hh"
#include "globals.hh"

class G4LatticePhysical;
class G4ParticleDefinition;
class G4Track;


namespace G4CMP {
  template <class T> G4int sign(T val) {
    return (T(0) < val) - (val < T(0));
  }

  // Select phonon mode randomly from density of states
  G4int ChoosePhononPolarization(const G4LatticePhysical* lattice);
  G4int ChoosePhononPolarization(G4double Ldos, G4double STdos, G4double FTdos);

  // Identify G4CMP particle categories
  G4bool IsPhonon(const G4Track* track);
  G4bool IsElectron(const G4Track* track);
  G4bool IsHole(const G4Track* track);
  G4bool IsChargeCarrier(const G4Track* track) {
    return (IsElectron(track) || IsHole(track));
  }

  G4bool IsPhonon(const G4ParticleDefinition* pd);
  G4bool IsElectron(const G4ParticleDefinition* pd);
  G4bool IsHole(const G4ParticleDefinition* pd);
  G4bool IsChargeCarrier(const G4ParticleDefinition* pd) {
    return (IsElectron(pd) || IsHole(pd));
  }
}

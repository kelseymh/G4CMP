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
//
// 20170602  Provide call-by-reference versions of track identity functions
// 20170802  Provide scale factor argument to ChooseWeight functions
// 20170928  Replace "polarization" with "mode"
// 20190906  Add function to get process associated with particle
// 20220816  Move RandomIndex function from SecondaryProduction
// 20220921  G4CMP-319 -- Add utilities for thermal (Maxwellian) distributions

#ifndef G4CMPUtils_hh
#define G4CMPUtils_hh 1

#include "G4ThreeVector.hh"
#include "globals.hh"

class G4CMPElectrodeHit;
class G4LatticePhysical;
class G4ParticleDefinition;
class G4Step;
class G4Track;
class G4VProcess;


namespace G4CMP {
  template <class T> G4int sign(T val) {
    return (T(0) < val) - (val < T(0));
  }

  // Identify G4CMP particle categories
  G4bool IsPhonon(const G4Track* track);
  G4bool IsElectron(const G4Track* track);
  G4bool IsHole(const G4Track* track);
  G4bool IsChargeCarrier(const G4Track* track);

  G4bool IsPhonon(const G4Track& track);
  G4bool IsElectron(const G4Track& track);
  G4bool IsHole(const G4Track& track);
  G4bool IsChargeCarrier(const G4Track& track);

  G4bool IsPhonon(const G4ParticleDefinition* pd);
  G4bool IsElectron(const G4ParticleDefinition* pd);
  G4bool IsHole(const G4ParticleDefinition* pd);
  G4bool IsChargeCarrier(const G4ParticleDefinition* pd);

  G4bool IsPhonon(const G4ParticleDefinition& pd);
  G4bool IsElectron(const G4ParticleDefinition& pd);
  G4bool IsHole(const G4ParticleDefinition& pd);
  G4bool IsChargeCarrier(const G4ParticleDefinition& pd);

  // Select phonon mode randomly from density of states
  G4int ChoosePhononPolarization(const G4LatticePhysical* lattice);
  G4int ChoosePhononPolarization(G4double Ldos, G4double STdos, G4double FTdos);

  // Randomly choose valley for charge carrier
  G4int ChooseValley(const G4LatticePhysical* lattice);

  // Throw biasing decision for particle production and return weight
  // NOTE:  biasScale < 0. means to use "primary generator" scaling
  G4double ChooseWeight(const G4ParticleDefinition* pd, G4double biasScale=-1.);
  G4double ChoosePhononWeight(G4double biasScale=-1.);
  G4double ChooseChargeWeight(G4double biasScale=-1.);

  // Create a Hit from a G4Step. Less error prone to use this helper.
  void FillHit(const G4Step*, G4CMPElectrodeHit*);

  // Phonons reflect difusively from surfaces.
  G4ThreeVector LambertReflection(const G4ThreeVector& surfNorm);

  // Test that a phonon's wave vector relates to an inward velocity.
  G4bool PhononVelocityIsInward(const G4LatticePhysical* lattice, G4int mode,
                                const G4ThreeVector& waveVector,
                                const G4ThreeVector& surfNorm);

  // Thermal distributions, useful for handling phonon thermalization
  G4double MaxwellBoltzmannPDF(G4double temperature, G4double energy);
  G4double ChooseThermalEnergy(G4double temperature);
  G4double ChooseThermalEnergy(const G4LatticePhysical* lattice);

  G4bool IsThermalized(G4double temperature, G4double energy);
  G4bool IsThermalized(G4double energy);	// Use G4CMPConfigManager temp.
  G4bool IsThermalized(const G4LatticePhysical* lattice, G4double energy);
  G4bool IsThermalized(const G4Track* track);
  inline G4bool IsThermalized(const G4Track& t) { return IsThermalized(&t); }

  // Search particle's processes for specified name
  G4VProcess* FindProcess(const G4ParticleDefinition* pd, const G4String& pname);

  // Generate integer random value [0, imax), used to shuffle vectors
  size_t RandomIndex(size_t imax);
}

#endif	/* G4CMPUtils_hh */

/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
// class G4CMPPhononElectrode, subclass of G4CMPVElectrodePattern
//
// Class description:
//
// This provides an interface to G4CMPKaplanQP, to simulate the response
// of a phonon energy-absorbing electrode made of a thin-film superconductor.
//
// The associated border surface must have a G4MaterialPropertiesTable
// registered, containing at least all the following parameters (see also
// README.md):
//
// | Property Key        | Definition                   | Example value (Al) |
// |---------------------|------------------------------|--------------------|
// | filmThickness       | Thickness of film            | 600.*nm            |
// | gapEnergy           | Bandgap of film material     | 173.715e-6*eV      |
// | lowQPLimit          | Minimum bandgap multiple     | 3.                 |
// | phononLifetime      | Phonon lifetime at 2*bandgap | 242.*ps            |
// | phononLifetimeSlope | Lifetime vs. energy          | 0.29               |
// | vSound              | Speed of sound in film       | 3.26*km/s          |
//
// In addition, for sensors containing a "direct absorber" (such as a TES),
// The property key "subgapAbsorption" may be set with the probability to
// directly absorb phonons below 2*bandgap.
// 
// 20221006  M. Kelsey -- Adapted from SuperCDMS simulation version

#ifndef G4CMPPhononElectrode_hh
#define G4CMPPhononElectrode_hh 1

#include "G4CMPVElectrodePattern.hh"

class G4CMPKaplanQP;
class G4ParticleChange;
class G4Step;
class G4Track;


class G4CMPPhononElectrode : public G4CMPVElectrodePattern {
public:
  G4CMPPhononElectrode();
  virtual ~G4CMPPhononElectrode();

  // Implement cloning function along with default copiers
  G4CMPPhononElectrode(const G4CMPPhononElectrode&) = default;
  G4CMPPhononElectrode(G4CMPPhononElectrode&&) = default;
  G4CMPPhononElectrode& operator=(const G4CMPPhononElectrode&) = default;
  G4CMPPhononElectrode& operator=(G4CMPPhononElectrode&&) = default;

  virtual G4CMPVElectrodePattern* Clone() const {
    return new G4CMPPhononElectrode(*this);
  }

  // Assumes that user has configured a border surface only at sensor pads
  virtual G4bool IsNearElectrode(const G4Step&) const;

  virtual void AbsorbAtElectrode(const G4Track&,
                                 const G4Step&,
                                 G4ParticleChange&) const;

protected:
  // Record energy deposition and re-emitted energies as secondary phonons
  void ProcessAbsorption(const G4Track& track, const G4Step& step,
			 G4double EDep, G4ParticleChange& particleChange) const;

  // Reflect phonon by changing direction randomly, no energy deposition
  void ProcessReflection(const G4Track& track, const G4Step& step,
			 G4ParticleChange& particleChange) const;

  // NOTE: "Mutable" because AbsorbAtElectrode() function is const
  mutable G4CMPKaplanQP* kaplanQP;	// Create instance of QET simulator
  mutable std::vector<G4double> phononEnergies;		// Reusable buffer
};

#endif

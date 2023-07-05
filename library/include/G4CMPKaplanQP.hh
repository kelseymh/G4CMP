/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPKaplanQP.hh
/// \brief Grouping of free standing functions that relate to the
/// creation and energy calculations of quasi-particle downconversion
/// by phonons breaking Cooper pairs in superconductors.
//
// $Id$
//
// 20200616  M. Kelsey -- Reimplement as class, keeping "KaplanPhononQP"
//		interface for migration.
// 20200618  G4CMP-212: Add optional parameter for below-bandgap phonons
//		to be absorbed in the superconducting film.
// 20200626  G4CMP-215: Add function to encapsulate below-bandgap absorption.
// 20200627  In *EnergyRand(), move PDF expressions to functions; eliminate
//		mutation of E argument in PhononEnergyRand().
// 20200701  G4CMP-217: New function to handle QP energy absorption below
//		minimum for QP -> phonon -> new QP pair chain (3*bandgap).
// 20201109  Add diagnostic text file (like downconversion and Luke).

#ifndef G4CMPKaplanQP_hh
#define G4CMPKaplanQP_hh 1

#include "G4Types.hh"
#include <fstream>
#include <vector>

class G4MaterialPropertiesTable;


// This is the main function for the Kaplan quasiparticle downconversion
// process. Based on the energy of the incoming phonon and the properties
// of the superconductor, we return the total energy deposited as well
// as fill a vector of energies that correspond to newly created phonons
// that are emitted back into the crystal.
namespace G4CMP {
  G4double KaplanPhononQP(G4double energy,
			  G4MaterialPropertiesTable* prop,
			  std::vector<G4double>& reflectedEnergies);
}

class G4CMPKaplanQP {
public:
  G4CMPKaplanQP(G4MaterialPropertiesTable* prop, G4int vb=0);
  virtual ~G4CMPKaplanQP();

  // Turn on diagnostic messages
  void SetVerboseLevel(G4int vb) { verboseLevel = vb; }
  G4int GetVerboseLevel() const { return verboseLevel; }

  // Configure thin film (QET, metalization, etc.) for phonon absorption
  void SetFilmProperties(G4MaterialPropertiesTable* prop);

  // Do absorption on sensor/metalization film
  // Returns absorbed energy, fills list of re-emitted phonons
  G4double AbsorbPhonon(G4double energy,
			std::vector<G4double>& reflectedEnergies) const;

protected:
  // Compute the probability of a phonon reentering the crystal without breaking
  // any Cooper pairs.
  G4double CalcEscapeProbability(G4double energy,
				 G4double thicknessFrac) const;

  // Model the phonons (phonEnergies) breaking Cooper pairs into quasiparticles
  // (qpEnergies).
  G4double CalcQPEnergies(std::vector<G4double>& phonEnergies,
			  std::vector<G4double>& qpEnergies) const;
  
  // Model the quasiparticles (qpEnergies) emitting phonons (phonEnergies) in
  // the superconductor.
  G4double CalcPhononEnergies(std::vector<G4double>& phonEnergies,
			      std::vector<G4double>& qpEnergies) const;
  
  // Calculate energies of phonon tracks that have reentered the crystal.
  void CalcReflectedPhononEnergies(std::vector<G4double>& phonEnergies,
				   std::vector<G4double>& reflectedEnergies) const;

  // Compute probability of absorbing phonon below Cooper-pair breaking
  // NOTE:  Function should ONLY be called for energy < 2.*gapEnergy
  G4double CalcSubgapAbsorption(G4double energy,
				std::vector<G4double>& keepEnergies) const;

  // Handle absorption of quasiparticle energies below Cooper-pair breaking
  // If qpEnergy < 3*Delta, radiate a phonon, absorb bandgap minimum
  G4double CalcQPAbsorption(G4double energy,
			    std::vector<G4double>& phonEnergies,
			    std::vector<G4double>& qpEnergies) const;

  // Compute quasiparticle energy distribution from broken Cooper pair.
  G4double QPEnergyRand(G4double Energy) const;
  G4double QPEnergyPDF(G4double E, G4double x) const;
  
  // Compute phonon energy distribution from quasiparticle in superconductor.
  G4double PhononEnergyRand(G4double Energy) const;
  G4double PhononEnergyPDF(G4double E, G4double x) const;

  // Encapsulate below-bandgap logic
  G4bool IsSubgap(G4double energy) const { return energy < 2.*gapEnergy; }

private:
  G4int verboseLevel;		// For diagnostic messages

  G4MaterialPropertiesTable* filmProperties;
  G4double filmThickness;	// Quantities extracted from properties table
  G4double gapEnergy;		// Bandgap energy (delta)
  G4double lowQPLimit;		// Minimum X*delta to keep as a quasiparticle
  G4double subgapAbsorption;	// Probability to absorb energy below bandgap
  G4double phononLifetime;	// Lifetime of phonons in film at 2*delta
  G4double phononLifetimeSlope;	// Energy dependence of phonon lifetime
  G4double vSound;		// Speed of sound in film

  mutable std::ofstream output;		// Diagnostic output under G4CMP_DEBUG
};

#endif	/* G4CMPKaplanQP_hh */

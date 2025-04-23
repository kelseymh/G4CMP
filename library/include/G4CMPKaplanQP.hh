/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPKaplanQP.hh
/// \brief Grouping of free standing functions that relate to the
/// creation and energy calculations of quasi-particle downconversion
/// by phonons breaking Cooper pairs in superconductors.
///
/// This code implements a "lumped" version of Kaplan's model for
/// quasiparticle-phonon interactions in superconducting films,
/// S.B.Kaplan et al., Phys.Rev.B14 (1976).
///
/// If the thin-film parameters are set from a MaterialPropertiesTable,
/// the table must contain the first five of the following entries:
///
/// | Property Key        | Definition                   | Example value (Al) |
/// |---------------------|------------------------------|--------------------|
/// | filmThickness       | Thickness of film            | 600.*nm            |
/// | vSound              | Speed of sound in film       | 3.26*km/s          |
/// | gapEnergy           | Bandgap of film material     | 173.715e-6*eV      |
/// | phononLifetime      | Phonon lifetime at 2*bandgap | 242.*ps            |
/// | phononLifetimeSlope | Lifetime vs. energy          | 0.29               |
/// |                     |                              |                    |
/// | lowQPLimit          | Minimum QP energy to radiate phonons | 3.         |
/// | highQPLimit         | Maximum energy to create QPs | 10.                |
/// | subgapAbsorption    | Absorption below 2*bandgap   | 0.03 (optional)    |
/// | absorberGap         | Bandgap of "subgap absorber" | 15e-6*eV (W)       |
/// | absorberEff         | QP absorption efficiency     | 0.3                |
/// | absorberEffSlope    | Efficiency vs. energy        | 0.                 |
/// | temperature         | Temperature of film          | 0.05e-3*K          |
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
// 20220928  G4CMP-323: Add bandgap of secondary absorber (quasiparticle trap)
//		Add direct-setting functions for configuration parameters,
//		and function to test whether parameters have been set.
// 20221006  G4CMP-330: Add temperature parameter with Set function.
// 20221102  G4CMP-314: Add energy dependent efficiency for QP absorption.
// 20221127  G4CMP-347: Add highQPLimit to split incident phonons
// 20221201  G4CMP-345: Rename "CalcSubgapAbs" to "CalcDirectAbs", split into
//		new DoDirectAbsorption() boolean test.
// 20240502  G4CMP-344: Reusable vector buffers to avoid memory churn.
// 20240502  G4CMP-379: Add Fermi-Dirac thermal probability for QP energies.

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

  // Do absorption on sensor/metalization film
  // Returns absorbed energy, fills list of re-emitted phonons
  G4double AbsorbPhonon(G4double energy,
			std::vector<G4double>& reflectedEnergies) const;

  // Set temperature for use by thermalization functions
  void SetTemperature(G4double temp) { temperature = temp; }

  // Configure thin film (QET, metalization, etc.) for phonon absorption
  void SetFilmProperties(G4MaterialPropertiesTable* prop);

  // Alternative configuration without properties table
  void SetFilmThickness(G4double value)       { filmThickness = value; }
  void SetGapEnergy(G4double value)           { gapEnergy = value; }
  void SetLowQPLimit(G4double value)          { lowQPLimit = value; }
  void SetHighQPLimit(G4double value)         { highQPLimit = value; }
  void SetSubgapAbsorption(G4double value)    { directAbsorption = value; }
  void SetDirectAbsorption(G4double value)    { directAbsorption = value; }
  void SetAbsorberGap(G4double value)         { absorberGap = value; }
  void SetAbsorberEff(G4double value)         { absorberEff = value; }
  void SetAbsorberEffSlope(G4double value)    { absorberEffSlope = value; }
  void SetPhononLifetime(G4double value)      { phononLifetime = value; }
  void SetPhononLifetimeSlope(G4double value) { phononLifetimeSlope = value; }
  void SetVSound(G4double value)              { vSound = value; }

protected:
  // Check that the five required parameters are set to meaningful values
  G4bool ParamsReady() const {
    return (filmThickness > 0. && gapEnergy >= 0. && vSound > 0. &&
	    phononLifetime > 0. && phononLifetimeSlope >= 0.);
  }

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

  // Compute probability of phonon collection directly on absorber (TES)
  G4bool DoDirectAbsorption(G4double energy) const;
  G4double CalcDirectAbsorption(G4double energy,
				std::vector<G4double>& keepEnergies) const;

  // Handle absorption of quasiparticle energies below Cooper-pair breaking
  // If qpEnergy < 3*Delta, radiate a phonon, absorb bandgap minimum
  G4double CalcQPAbsorption(G4double energy,
			    std::vector<G4double>& phonEnergies,
			    std::vector<G4double>& qpEnergies) const;
			    
  // Handle quasiparticle energy-dependent absorption efficiency
  G4double CalcQPEfficiency(G4double qpE) const;

  // Compute quasiparticle energy distribution from broken Cooper pair.
  G4double QPEnergyRand(G4double Energy) const;
  G4double QPEnergyPDF(G4double E, G4double x) const;
  G4double ThermalPDF(G4double E) const;

  // Compute phonon energy distribution from quasiparticle in superconductor.
  G4double PhononEnergyRand(G4double Energy) const;
  G4double PhononEnergyPDF(G4double E, G4double x) const;

  // Encapsulate below-bandgap logic
  G4bool IsSubgap(G4double energy) const { return (energy < 2.*gapEnergy); }
  G4bool DirectAbsorb(G4double energy) const {
    return (IsSubgap(energy) && energy > 2.*absorberGap);
  }

  // Write summary of interaction to output "kaplanqp_stats" file
  void ReportAbsorption(G4double energy, G4double EDep,
			const std::vector<G4double>& reflectedEnergies) const;

private:
  G4int verboseLevel;			// For diagnostic messages
  mutable G4bool keepAllPhonons;	// Copy of flag KeepKaplanPhonons()

  G4MaterialPropertiesTable* filmProperties;
  G4double filmThickness;	// Quantities extracted from properties table
  G4double gapEnergy;		// Bandgap energy (delta)
  G4double lowQPLimit;		// Minimum X*delta to keep as a quasiparticle
  G4double highQPLimit;		// Maximum X*delta to create QP from phonon
  G4double directAbsorption;	// Probability to collect energy directly (TES)
  G4double absorberGap;		// Bandgap of secondary absorber material
  G4double absorberEff;         // Quasiparticle absorption efficiency
  G4double absorberEffSlope;    // Energy dependence of qp absorption efficiency
  G4double phononLifetime;	// Lifetime of phonons in film at 2*delta
  G4double phononLifetimeSlope;	// Energy dependence of phonon lifetime
  G4double vSound;		// Speed of sound in film
  G4double temperature;		// Ambient temperature of film (from lattice)

  // Temporary buffers for use within processing functions
  mutable std::vector<G4double> qpEnergyList;	// Active ("final") population
  mutable std::vector<G4double> phononEnergyList;
  mutable std::vector<G4double> newQPEnergies;	// Intermediate processing
  mutable std::vector<G4double> newPhonEnergies;

  mutable std::ofstream output;		// Diagnostic output under G4CMP_DEBUG
};

#endif	/* G4CMPKaplanQP_hh */

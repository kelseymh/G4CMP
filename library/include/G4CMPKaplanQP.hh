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
// 20250101  G4CMP-439: Create separate debugging file per worker thread;
//    add EventID and TrackID columns to debugging output.
// 20260710  G4CMP-647 -- Derive from new base class G4CMPVKaplanQP

#ifndef G4CMPKaplanQP_hh
#define G4CMPKaplanQP_hh 1

#include "G4CMPVKaplanQP.hh"
#include "G4Types.hh"
#include <fstream>
#include <vector>

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

class G4CMPKaplanQP : public G4CMPVKaplanQP {
public:
  G4CMPKaplanQP(G4MaterialPropertiesTable* prop, G4int vb=0);
  virtual ~G4CMPKaplanQP() {;}

  // Do absorption on sensor/metalization film
  // Returns absorbed energy, fills list of re-emitted phonons
  virtual G4double AbsorbPhonon(G4double energy,
			std::vector<G4double>& reflectedEnergies) const override;

protected:
  // Compute the probability of a phonon reentering the crystal without breaking
  // any Cooper pairs.
  virtual G4double CalcEscapeProbability(G4double energy,
				 G4double thicknessFrac) const override;

  // Model the phonons (phonEnergies) breaking Cooper pairs into quasiparticles
  // (qpEnergies).
  virtual G4double CalcQPEnergies(std::vector<G4double>& phonEnergies,
			  std::vector<G4double>& qpEnergies) const override;
  
  // Model the quasiparticles (qpEnergies) emitting phonons (phonEnergies) in
  // the superconductor.
  virtual G4double CalcPhononEnergies(std::vector<G4double>& phonEnergies,
			      std::vector<G4double>& qpEnergies) const override;
  
  // Calculate energies of phonon tracks that have reentered the crystal.
  virtual void CalcReflectedPhononEnergies(std::vector<G4double>& phonEnergies,
				   std::vector<G4double>& reflectedEnergies) const override;

  // Compute probability of phonon collection directly on absorber (TES)
  virtual G4bool DoDirectAbsorption(G4double energy) const override;
  virtual G4double CalcDirectAbsorption(G4double energy,
				std::vector<G4double>& keepEnergies) const override;

  // Handle absorption of quasiparticle energies below Cooper-pair breaking
  // If qpEnergy < 3*Delta, radiate a phonon, absorb bandgap minimum
  virtual G4double CalcQPAbsorption(G4double energy,
			    std::vector<G4double>& phonEnergies,
			    std::vector<G4double>& qpEnergies) const override;
			    
  // Handle quasiparticle energy-dependent absorption efficiency
  virtual G4double CalcQPEfficiency(G4double qpE) const override;
}; // class G4CMPKaplanQP

#endif	/* G4CMPKaplanQP_hh */

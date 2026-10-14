/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPVKaplanQP.hh
/// \brief Base class for phonon-qp interaction physics
///
/// This code implements a base class to control the properties of
/// phonon-qp interactions. Most of the code is extracted from the original
/// G4CMPKaplanQP implementation.
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
// 20260710  G4CMP-647 -- Implement virtual base class for phonon-qp
// interactions

#ifndef G4CMPVKaplanQP_hh
#define G4CMPVKaplanQP_hh 1

#include "G4Types.hh"
#include <fstream>
#include <vector>

class G4MaterialPropertiesTable;

class G4CMPVKaplanQP {
public:
  G4CMPVKaplanQP(G4MaterialPropertiesTable* prop, G4int vb=0);
  virtual ~G4CMPVKaplanQP();

  // Turn on diagnostic messages
  virtual void SetVerboseLevel(G4int vb) { verboseLevel = vb; }
  G4int GetVerboseLevel() const { return verboseLevel; }

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

  // Do absorption on sensor/metalization film
  // Returns absorbed energy, fills list of re-emitted phonons
  virtual G4double AbsorbPhonon(G4double energy,
      std::vector<G4double>& reflectedEnergies) const = 0;

protected:
  // Check that the five required parameters are set to meaningful values
  virtual G4bool ParamsReady() const {
    return (filmThickness > 0. && gapEnergy > 0. && vSound > 0. &&
      phononLifetime > 0. && phononLifetimeSlope >= 0.);
  }

  // Compute the probability of a phonon reentering the crystal without breaking
  // any Cooper pairs.
  virtual G4double CalcEscapeProbability(G4double energy,
         G4double thicknessFrac) const = 0;

  // Model the phonons (phonEnergies) breaking Cooper pairs into quasiparticles
  // (qpEnergies).
  virtual G4double CalcQPEnergies(std::vector<G4double>& phonEnergies,
        std::vector<G4double>& qpEnergies) const = 0;
  
  // Model the quasiparticles (qpEnergies) emitting phonons (phonEnergies) in
  // the superconductor.
  virtual G4double CalcPhononEnergies(std::vector<G4double>& phonEnergies,
            std::vector<G4double>& qpEnergies) const = 0;
  
  // Calculate energies of phonon tracks that have reentered the crystal.
  virtual void CalcReflectedPhononEnergies(std::vector<G4double>& phonEnergies,
           std::vector<G4double>& reflectedEnergies) const = 0;

  // Compute probability of phonon collection directly on absorber (TES)
  virtual G4bool DoDirectAbsorption(G4double energy) const = 0;
  virtual G4double CalcDirectAbsorption(G4double energy,
        std::vector<G4double>& keepEnergies) const = 0;

  // Handle absorption of quasiparticle energies below Cooper-pair breaking
  // If qpEnergy < 3*Delta, radiate a phonon, absorb bandgap minimum
  virtual G4double CalcQPAbsorption(G4double energy,
          std::vector<G4double>& phonEnergies,
          std::vector<G4double>& qpEnergies) const = 0;
          
  // Handle quasiparticle energy-dependent absorption efficiency
  virtual G4double CalcQPEfficiency(G4double qpE) const = 0;

  // Compute quasiparticle energy distribution from broken Cooper pair.
  G4double QPEnergyRand(G4double Energy) const;
  G4double QPEnergyPDF(G4double E, G4double x) const;

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

protected:
  G4int verboseLevel;     // For diagnostic messages
  mutable G4bool keepAllPhonons;  // Copy of flag KeepKaplanPhonons()

  G4MaterialPropertiesTable* filmProperties;
  G4double filmThickness; // Quantities extracted from properties table
  G4double gapEnergy;   // Bandgap energy (delta)
  G4double lowQPLimit;    // Minimum X*delta to keep as a quasiparticle
  G4double highQPLimit;   // Maximum X*delta to create QP from phonon
  G4double directAbsorption;  // Probability to collect energy directly (TES)
  G4double absorberGap;   // Bandgap of secondary absorber material
  G4double absorberEff;         // Quasiparticle absorption efficiency
  G4double absorberEffSlope;    // Energy dependence of qp absorption efficiency
  G4double phononLifetime;  // Lifetime of phonons in film at 2*delta
  G4double phononLifetimeSlope; // Energy dependence of phonon lifetime
  G4double vSound;    // Speed of sound in film
  G4double temperature;   // Ambient temperature of film (from lattice)

  // Temporary buffers for use within processing functions
  mutable std::vector<G4double> qpEnergyList; // Active ("final") population
  mutable std::vector<G4double> phononEnergyList;
  mutable std::vector<G4double> newQPEnergies;  // Intermediate processing
  mutable std::vector<G4double> newPhonEnergies;

  mutable std::ofstream output;   // Diagnostic output under G4CMP_DEBUG
}; // class G4CMPVKaplanQP

#endif  /* G4CMPVKaplanQP_hh */

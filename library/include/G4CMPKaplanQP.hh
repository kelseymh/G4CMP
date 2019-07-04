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

#ifndef G4CMPKaplanQP_hh
#define G4CMPKaplanQP_hh 1

#include "globals.hh"
#include <vector>

class G4MaterialPropertiesTable;

namespace G4CMP {

// This is the main function for the Kaplan quasiparticle downconversion
// process. Based on the energy of the incoming phonon and the properties
// of the superconductor, we return the total energy deposited as well
// as fill a vector of energies that correspond to newly created phonons
// that are emitted back into the crystal.
G4double KaplanPhononQP(G4double energy,
                        G4MaterialPropertiesTable* prop,
                        std::vector<G4double>& reflectedEnergies);

// Compute the probability of a phonon reentering the crystal without breaking
// any Cooper pairs.
G4double CalcEscapeProbability(G4double energy,
                               G4double thicknessFrac,
                               G4MaterialPropertiesTable* prop);

// Model the phonons (phonEnergies) breaking Cooper pairs into quasiparticles
// (qpEnergies).
G4double CalcQPEnergies(G4double gapEnergy,
                        G4double lowQPLimit,
                        std::vector<G4double>& phonEnergies,
                        std::vector<G4double>& qpEnergies);

// Model the quasiparticles (qpEnergies) emitting phonons (phonEnergies) in
// the superconductor.
G4double CalcPhononEnergies(G4double gapEnergy,
                            G4double lowQPLimit,
                            std::vector<G4double>& phonEnergies,
                            std::vector<G4double>& qpEnergies);

// Calculate energies of phonon tracks that have reentered the crystal.
void CalcReflectedPhononEnergies(G4MaterialPropertiesTable* prop,
                                 std::vector<G4double>& phonEnergies,
                                 std::vector<G4double>& reflectedEnergies);

// Compute quasiparticle energy distribution from broken Cooper pair.
G4double QPEnergyRand(G4double gapEnergy, G4double Energy);

// Compute phonon energy distribution from quasiparticle in superconductor.
G4double PhononEnergyRand(G4double gapEnergy, G4double& Energy);
}

#endif

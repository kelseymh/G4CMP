/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPKaplanQP.cc
/// \brief Grouping of free standing functions that relate to the
/// creation and energy calculations of quasi-particle downconversion
/// by phonons breaking Cooper pairs in superconductors.
//
// $Id$

#include "G4CMPKaplanQP.hh"
#include "G4MaterialPropertiesTable.hh"
#include "Randomize.hh"

G4double G4CMP::KaplanPhononQP(G4double energy,
                               G4MaterialPropertiesTable* prop,
                               std::vector<G4double>& reflectedEnergies) {
  if (reflectedEnergies.size() > 0) {
    G4Exception("G4CMP::KaplanPhononQP()", "G4CMP007",
                JustWarning, "Passed a nonempty vector.");
  }

  // Check that the MaterialPropertiesTable has everything we need. If it came
  // from a G4CMPSurfaceProperty, then it will be fine.
  if (!(prop->ConstPropertyExists("gapEnergy") &&
        prop->ConstPropertyExists("lowQPLimit") &&
        prop->ConstPropertyExists("phononLifetime") &&
        prop->ConstPropertyExists("phononLifetimeSlope") &&
        prop->ConstPropertyExists("vSound") &&
        prop->ConstPropertyExists("filmThickness"))) {
    G4Exception("G4CMP::KaplanPhononQP()", "G4CMP001",
                RunMustBeAborted,
                "Insufficient info in MaterialPropertiesTable.");
  }

  G4double gapEnergy     = prop->GetConstProperty("gapEnergy");
  G4double lowQPLimit    = prop->GetConstProperty("lowQPLimit");

  // For the phonon to not break a Cooper pair, it must go 2*thickness,
  // assuming it goes exactly along the thickness direction, which is an
  // approximation.
  G4double frac = 2.0;
  G4double phononEscapeProb = CalcEscapeProbability(energy, frac, prop);

  G4double EDep = 0.;
  if (energy > 2.0*gapEnergy && G4UniformRand() > phononEscapeProb) {
    std::vector<G4double> qpEnergies;
    std::vector<G4double> phonEnergies{energy};
    while (qpEnergies.size() > 0 || phonEnergies.size() > 0) {
      if (phonEnergies.size() > 0) {
        // Partition the phonons' energies into quasi-particles according to
        // a PDF defined in CalcQPEnergies().
        // NOTE: Both energy vectors mutate.
        EDep += CalcQPEnergies(gapEnergy, lowQPLimit, phonEnergies, qpEnergies);
      }
      if (qpEnergies.size() > 0) {
        // Quasiparticles can also excite phonons.
        // NOTE: Both energy vectors mutate.
        EDep += CalcPhononEnergies(gapEnergy, lowQPLimit, phonEnergies, qpEnergies);
      }
      if (phonEnergies.size() > 0) {
        // Some phonons will escape back into the crystal.
        // NOTE: Both energy vectors mutate.
        CalcReflectedPhononEnergies(prop, phonEnergies, reflectedEnergies);
      }
    }
  } else {
    reflectedEnergies.push_back(energy);
  }

  return EDep;
}

G4double G4CMP::CalcEscapeProbability(G4double energy,
                                      G4double thicknessFrac,
                                      G4MaterialPropertiesTable* prop) {
  // Grab all necessary properties from the material.
  G4double gapEnergy = prop->GetConstProperty("gapEnergy");
  G4double phononLifetime = prop->GetConstProperty("phononLifetime");
  G4double phononLifetimeSlope = prop->GetConstProperty("phononLifetimeSlope");
  G4double vSound = prop->GetConstProperty("vSound");
  G4double thickness = prop->GetConstProperty("filmThickness");

  G4double mfp = vSound * phononLifetime /
                 (1. + phononLifetimeSlope * (energy/gapEnergy - 2.));
  return std::exp(-2.* thicknessFrac * thickness/mfp);
}

G4double G4CMP::CalcQPEnergies(G4double gapEnergy,
                               G4double lowQPLimit,
                               std::vector<G4double>& phonEnergies,
                               std::vector<G4double>& qpEnergies) {
  // Each phonon gives all of its energy to the qp pair it breaks.
  G4double EDep = 0.;
  for (G4double E: phonEnergies) {
    G4double qpE = QPEnergyRand(gapEnergy, E);
    if (qpE >= lowQPLimit*gapEnergy) {
      qpEnergies.push_back(qpE);
    } else {
      EDep += qpE;
    }

    if (E-qpE >= lowQPLimit*gapEnergy) {
      qpEnergies.push_back(E-qpE);
    } else {
      EDep += E-qpE;
    }
  }

  phonEnergies.clear();
  return EDep;
}


G4double G4CMP::CalcPhononEnergies(G4double gapEnergy,
                                   G4double lowQPLimit,
                                   std::vector<G4double>& phonEnergies,
                                   std::vector<G4double>& qpEnergies) {
  // NOTE: Phonons with low energy will not be seen by the detector, so we
  // don't record those energies and just "lose" those phonons.
  // Have a reference in for loop b/c qp doesn't give all of its energy away.
  G4double EDep = 0.;
  std::vector<G4double> newQPEnergies;
  for (G4double& E: qpEnergies) {
    // NOTE: E mutates in PhononEnergyRand.
    G4double phonE = PhononEnergyRand(gapEnergy, E);
    if (phonE >= 2.0*gapEnergy) {
      phonEnergies.push_back(phonE);
    }
    if (E >= lowQPLimit*gapEnergy) {
      newQPEnergies.push_back(E);
    } else {
      EDep += E;
    }
  }

  qpEnergies.swap(newQPEnergies);
  return EDep;
}

void G4CMP::CalcReflectedPhononEnergies(G4MaterialPropertiesTable* prop,
                                        std::vector<G4double>& phonEnergies,
                                        std::vector<G4double>& reflectedEnergies) {
  // There is a 50% chance that a phonon is headed away from (toward) substrate
  std::vector<G4double> newPhonEnergies;
  for (G4double E : phonEnergies) {
    // frac = 1.5 for phonons headed away from the subst. 0.5 for toward.
    // This assumes that, on average, the phonons are spawned at the center
    // of the superconductor, which is likely not true.
    G4double frac = (G4UniformRand() < 0.5 ? 0.5 : 1.5);
    if (G4UniformRand() < CalcEscapeProbability(E, frac, prop)) {
      reflectedEnergies.push_back(E);
    } else {
      newPhonEnergies.push_back(E);
    }
  }
  phonEnergies.swap(newPhonEnergies);
}

G4double G4CMP::QPEnergyRand(G4double gapEnergy, G4double Energy) {
  // PDF is not integrable, so we can't do an inverse transform sampling.
  // Instead, we'll do a rejection method.
  //
  // PDF(E') = (E'*(Energy - E') + gapEnergy*gapEnergy)
  //           /
  //           sqrt((E'*E' - gapEnergy*gapEnergy) *
  //                ((Energy - E')*(Energy - E') - gapEnergy*gapEnergy));
  // The shape of the PDF is like a U, so the max values are at the endpoints:
  // E' = gapEnergy and E' = Energy - gapEnergy

  // Add buffer so first/last bins don't give zero denominator in pdfSum
  const G4double BUFF = 1000.;
  G4double xmin = gapEnergy + (Energy-2.*gapEnergy)/BUFF;
  G4double xmax = gapEnergy + (Energy-2.*gapEnergy)*(BUFF-1.)/BUFF;

  G4double ymax = (xmin*(Energy-xmin) + gapEnergy*gapEnergy)
                  /
                  sqrt((xmin*xmin - gapEnergy*gapEnergy) *
                       ((Energy-xmin)*(Energy-xmin) - gapEnergy*gapEnergy));

  G4double ytest = G4UniformRand()*ymax;
  G4double xtest = G4UniformRand()*(xmax-xmin) + xmin;
  while (ytest > (xtest*(Energy-xtest) + gapEnergy*gapEnergy)
                  /
                  sqrt((xtest*xtest - gapEnergy*gapEnergy) *
                       ((Energy-xtest)*(Energy-xtest) - gapEnergy*gapEnergy))) {
    ytest = G4UniformRand()*ymax;
    xtest = G4UniformRand()*(xmax-xmin) + xmin;
  }

  return xtest;
/*
  // PDF is not integrable, so we can't do an inverse transform sampling.
  // PDF is shaped like a capital U, so a rejection method would be slow.
  // Let's numerically calculate a CDF and throw a random to find E.

  const G4int BINS = 1000;
  G4double energyArr[BINS];
  G4double cdf[BINS];
  G4double pdfSum = 0.;
  for (size_t i = 0; i < BINS; ++i) {
    // Add 1 to i so first bin doesn't give zero denominator in pdfSum
    // Add 1 to BINS so last bin doesn't give zero denominator in pdfSum
    energyArr[i] = gapEnergy + (Energy-2.*gapEnergy) * (i+1.)/(BINS+1.);
    pdfSum += (energyArr[i]*(Energy - energyArr[i]) + gapEnergy*gapEnergy)
              /
              sqrt((energyArr[i]*energyArr[i] - gapEnergy*gapEnergy) *
                   ((Energy - energyArr[i])*(Energy - energyArr[i]) -
                    gapEnergy*gapEnergy));
    cdf[i] = pdfSum;
  }

  G4double u = G4UniformRand();

  size_t index = 0;
  for (; index < BINS; ++index) { //Combine normalization and search loops
    if (cdf[index]/pdfSum >= u)
      break;
  }

  return energyArr[index];
*/
}

G4double G4CMP::PhononEnergyRand(G4double gapEnergy, G4double& Energy) {
  // PDF is not integrable, so we can't do an inverse transform sampling.
  // Instead, we'll do a rejection method.
  //
  // PDF(E') = (E'*(Energy-E')*(Energy-E') * (E'-gapEnergy*gapEnergy/Energy))
  //           /
  //           sqrt((E'*E' - gapEnergy*gapEnergy);

  // Add buffer so first bin doesn't give zero denominator in pdfSum
  const G4double BUFF = 1000.;
  G4double xmin = gapEnergy + gapEnergy/BUFF;
  G4double xmax = Energy;

  G4double ymax = (xmin*(Energy-xmin)*(Energy-xmin) *
                    (xmin-gapEnergy*gapEnergy/Energy)) /
                  sqrt(xmin*xmin - gapEnergy*gapEnergy);

  G4double ytest = G4UniformRand()*ymax;
  G4double xtest = G4UniformRand()*(xmax-xmin) + xmin;
  while (ytest > (xtest*(Energy-xtest)*(Energy-xtest) *
                    (xtest-gapEnergy*gapEnergy/Energy)) /
                  sqrt(xtest*xtest - gapEnergy*gapEnergy)) {
    ytest = G4UniformRand()*ymax;
    xtest = G4UniformRand()*(xmax-xmin) + xmin;
  }

  G4double phononE = Energy - xtest;
  Energy = xtest;
  return phononE;
}

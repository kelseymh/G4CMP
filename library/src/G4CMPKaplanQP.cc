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
//
// 20200616  M. Kelsey -- Reimplement as class, keeping "KaplanPhononQP"
//		interface for migration.
// 20200618  G4CMP-212: Add optional parameter for below-bandgap phonons
//		to be absorbed in the superconducting film.

#include "globals.hh"
#include "G4CMPKaplanQP.hh"
#include "G4CMPConfigManager.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <numeric>


// Global function for Kaplan quasiparticle downconversion.  Retained here
// temporarily for migration to factory class

G4double G4CMP::KaplanPhononQP(G4double energy,
                               G4MaterialPropertiesTable* prop,
                               std::vector<G4double>& reflectedEnergies) {
  // Reusable factory class, which can be reconfigured on each call
  static G4ThreadLocal G4CMPKaplanQP* theKaplanQP = new G4CMPKaplanQP(prop);
  theKaplanQP->SetFilmProperties(prop);

  // TEMPORARY:  Set verbosity using the global value
  theKaplanQP->SetVerboseLevel(G4CMPConfigManager::GetVerboseLevel());

  return theKaplanQP->AbsorbPhonon(energy, reflectedEnergies);
}


// Class constructor

G4CMPKaplanQP::G4CMPKaplanQP(G4MaterialPropertiesTable* prop, G4int vb)
  : verboseLevel(vb), filmProperties(0), filmThickness(0.), gapEnergy(0.),
    lowQPLimit(0.), subgapAbsorption(0.), phononLifetime(0.),
    phononLifetimeSlope(0.), vSound(0.) {
  SetFilmProperties(prop);
}

// Configure thin film (QET, metalization, etc.) for phonon absorption

void G4CMPKaplanQP::SetFilmProperties(G4MaterialPropertiesTable* prop) {
  if (!prop) {
    G4Exception("G4CMPKaplanQP::SetFilmProperties()", "G4CMP001",
                RunMustBeAborted, "Null MaterialPropertiesTable vector.");
  }

  // Check that the MaterialPropertiesTable has everything we need. If it came
  // from a G4CMPSurfaceProperty, then it will be fine.
  if (!(prop->ConstPropertyExists("gapEnergy") &&
        prop->ConstPropertyExists("lowQPLimit") &&
        prop->ConstPropertyExists("phononLifetime") &&
        prop->ConstPropertyExists("phononLifetimeSlope") &&
        prop->ConstPropertyExists("vSound") &&
        prop->ConstPropertyExists("filmThickness"))) {
    G4Exception("G4CMPKaplanQP::SetFilmProperties()", "G4CMP002",
		RunMustBeAborted,
                "Insufficient info in MaterialPropertiesTable.");
  }

  // Extract values from table here for convenience in functions
  if (filmProperties != prop) {
    filmThickness =       prop->GetConstProperty("filmThickness");
    gapEnergy =           prop->GetConstProperty("gapEnergy");
    lowQPLimit =          prop->GetConstProperty("lowQPLimit");
    phononLifetime =      prop->GetConstProperty("phononLifetime");
    phononLifetimeSlope = prop->GetConstProperty("phononLifetimeSlope");
    vSound =              prop->GetConstProperty("vSound");

    subgapAbsorption = (prop->ConstPropertyExists("subgapAbsorption")
			? prop->GetConstProperty("subgapAbsorption") : 0.);

    filmProperties = prop;
  }
}


// This is the main function for the Kaplan quasiparticle downconversion
// process. Based on the energy of the incoming phonon and the properties
// of the superconductor, we return the total energy deposited as well
// as fill a vector of energies that correspond to newly created phonons
// that are emitted back into the crystal.

G4double G4CMPKaplanQP::
AbsorbPhonon(G4double energy, std::vector<G4double>& reflectedEnergies) const {
  if (!filmProperties) {
    G4Exception("G4CMPKaplanQP::AbsorbPhonon()", "G4CMP001",
                RunMustBeAborted, "Null MaterialPropertiesTable vector.");
  }

  if (verboseLevel)
    G4cout << "G4CMPKaplanQP::AbsorbPhonon " << energy << G4endl;

  if (reflectedEnergies.size() > 0) {
    G4Exception("G4CMPKaplanQP::AbsorbPhonon", "G4CMP007", JustWarning,
                "Passed a nonempty reflectedEnergies vector.");
    // FIXME: Should we discard previous contents?
  }

  // For the phonon to not break a Cooper pair, it must go 2*thickness,
  // assuming it goes exactly along the thickness direction, which is an
  // approximation.
  G4double frac = 2.0;
  G4double phononEscapeProb = CalcEscapeProbability(energy, frac);

  // If phonon is not absorbed, reflect it back with no deposition
  if (energy <= 2.0*gapEnergy && G4UniformRand() <= phononEscapeProb) {
    if (verboseLevel>1) G4cout << " Not absorbed." << G4endl;

    reflectedEnergies.push_back(energy);
    return 0.;
  }

  // Phonon goes into superconductor and gets partitioned into
  // quasiparticles, new phonons, and absorbed energy
  G4double EDep = 0.;

  std::vector<G4double> qpEnergies;
  std::vector<G4double> phonEnergies{energy};
  while (qpEnergies.size() > 0 || phonEnergies.size() > 0) {
    if (phonEnergies.size() > 0) {
      // Partition the phonons' energies into quasi-particles according to
      // a PDF defined in CalcQPEnergies().
      // NOTE: Both energy vectors mutate.
      EDep += CalcQPEnergies(phonEnergies, qpEnergies);
    }
    if (qpEnergies.size() > 0) {
      // Quasiparticles can also excite phonons.
      // NOTE: Both energy vectors mutate.
      EDep += CalcPhononEnergies(phonEnergies, qpEnergies);
    }
    if (phonEnergies.size() > 0) {
      // Some phonons will escape back into the crystal.
      // NOTE: Both energy vectors mutate.
      CalcReflectedPhononEnergies(phonEnergies, reflectedEnergies);
    }
  }

  // Sanity check -- Reflected + Absorbed should equal input
  G4double ERefl = std::accumulate(reflectedEnergies.begin(),
				   reflectedEnergies.end(), 0.);
  if (verboseLevel>1)
    G4cout << " Reflected " << ERefl << "\n Absorbed " << EDep << G4endl;

  if (fabs(energy-ERefl-EDep)/energy > 1e-3) {
    G4cerr << "WARNING G4CMPKaplanQP lost " << (energy-ERefl+EDep)/eV << " eV"
	   << G4endl;
  }

  return EDep;
}

G4double G4CMPKaplanQP::CalcEscapeProbability(G4double energy,
					      G4double thicknessFrac) const {
  if (verboseLevel>1) {
    G4cout << "G4CMPKaplanQP::CalcEscapeProbability E " << energy
	   << " thickFrac " << thicknessFrac << G4endl;
  }

  // Compute energy-dependent mean free path for phonons in film
  if (gapEnergy <= 0.) return 1.;

  G4double mfp = vSound * phononLifetime /
                 (1. + phononLifetimeSlope * (energy/gapEnergy - 2.));

  if (verboseLevel>2) {
    G4cout << " mfp " << mfp << " returning "
	   << std::exp(-2.* thicknessFrac * filmThickness/mfp) << G4endl;
  }

  return std::exp(-2.* thicknessFrac * filmThickness/mfp);
}

G4double 
G4CMPKaplanQP::CalcQPEnergies(std::vector<G4double>& phonEnergies,
			      std::vector<G4double>& qpEnergies) const {
  if (verboseLevel>1) {
    G4cout << "G4CMPKaplanQP::CalcQPEnergies QPcut " << lowQPLimit*gapEnergy
	   << G4endl;
  }

  // Each phonon gives all of its energy to the qp pair it breaks.
  G4double EDep = 0.;
  for (G4double E: phonEnergies) {
    G4double qpE = QPEnergyRand(E);
    if (verboseLevel>2) G4cout << " phononE " << E << " qpE " << qpE << G4endl;

    if (qpE >= lowQPLimit*gapEnergy) {
      if (verboseLevel>2) G4cout << " Storing qpE in qpEnergies" << G4endl;
      qpEnergies.push_back(qpE);
    } else {
      EDep += qpE;
    }

    if (E-qpE >= lowQPLimit*gapEnergy) {
      if (verboseLevel>2) G4cout << " Storing E-qpE in qpEnergies" << G4endl;
      qpEnergies.push_back(E-qpE);
    } else {
      EDep += E-qpE;
    }
  }

  phonEnergies.clear();		// All phonons have been processed

  if (verboseLevel>1) G4cout << " returning EDep " << EDep << G4endl;
  return EDep;
}


G4double 
G4CMPKaplanQP::CalcPhononEnergies(std::vector<G4double>& phonEnergies,
				  std::vector<G4double>& qpEnergies) const {
  if (verboseLevel>1) {
    G4cout << "G4CMPKaplanQP::CalcPhononEnergies 2*gap " << 2.*gapEnergy
	   << " QPcut " << lowQPLimit*gapEnergy << G4endl;
  }

  // Have a reference in for loop b/c qp doesn't give all of its energy away.
  G4double EDep = 0.;
  std::vector<G4double> newQPEnergies;
  for (G4double& E: qpEnergies) {
    if (verboseLevel>2) G4cout << " qpE " << E;		// Report before change

    // NOTE: E mutates in PhononEnergyRand.
    G4double phonE = PhononEnergyRand(E);
    if (verboseLevel>2) G4cout << " phononE " << phonE << " E " << E << G4endl;

    if (phonE >= 2.0*gapEnergy) {
      if (verboseLevel>2) G4cout << " Store phonE in phonEnergies" << G4endl;
      phonEnergies.push_back(phonE);
    } else if (G4UniformRand() < subgapAbsorption) {
      if (verboseLevel>2) G4cout << " Deposit phonE as heat" << G4endl;
      EDep += phonE;
    } else {
      if (verboseLevel>2) G4cout << " Return phonE for reflection" << G4endl;
      phonEnergies.push_back(phonE);
    }

    if (E >= lowQPLimit*gapEnergy) {
      if (verboseLevel>2) G4cout << " Store E in qpEnergies" << G4endl;
      newQPEnergies.push_back(E);
    } else {
      EDep += E;
    }
  }

  if (verboseLevel>1)
    G4cout << " replacing qpEnergies, returning EDep " << EDep << G4endl;

  qpEnergies.swap(newQPEnergies);
  return EDep;
}

void G4CMPKaplanQP::
CalcReflectedPhononEnergies(std::vector<G4double>& phonEnergies,
                            std::vector<G4double>& reflectedEnergies) const {
  if (verboseLevel>1)
    G4cout << "G4CMPKaplanQP::CalcReflectedPhononEnergies " << G4endl;

  // There is a 50% chance that a phonon is headed away from (toward) substrate
  std::vector<G4double> newPhonEnergies;
  for (G4double E : phonEnergies) {
    if (verboseLevel>2) G4cout << " phononE " << E << G4endl;

    // frac = 1.5 for phonons headed away from the subst. 0.5 for toward.
    // This assumes that, on average, the phonons are spawned at the center
    // of the superconductor, which is likely not true.
    G4double frac = (G4UniformRand() < 0.5 ? 0.5 : 1.5);
    if (G4UniformRand() < CalcEscapeProbability(E, frac)) {
      if (verboseLevel>2) G4cout << " phononE got reflected" << G4endl;
      reflectedEnergies.push_back(E);
    } else {
      if (verboseLevel>2) G4cout << " phononE stays in film" << G4endl;
      newPhonEnergies.push_back(E);
    }
  }
  phonEnergies.swap(newPhonEnergies);
}

G4double G4CMPKaplanQP::QPEnergyRand(G4double Energy) const {
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

// NOTE:  Input "Energy" value is replaced with E-phononE on return
G4double G4CMPKaplanQP::PhononEnergyRand(G4double& Energy) const {
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

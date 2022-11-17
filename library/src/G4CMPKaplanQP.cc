/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPKaplanQP.cc
/// \brief Grouping of free standing functions that relate to the
/// creation and energy calculations of quasi-particle downconversion
/// by phonons breaking Cooper pairs in superconductors.
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
/// | lowQPLimit          | Minimum bandgap multiple     | 3.                 |
/// | subgapAbsorption    | Absorption below 2*bandgap   | 0.03 (optional)    |
/// | absorberGap         | Bandgap of "subgap absorber" | 15e-6*eV (W)       |
/// | absorberEff	  | QP absorption efficiency     | 0.3   	      |
/// | absorberEffSlope    | Efficiency vs. energy        | 0.                 |
/// | temperature         | Temperature of film          | 0.05e-3*K          |
//
// $Id$
//
// 20200616  M. Kelsey -- Reimplement as class, keeping "KaplanPhononQP"
//		interface for migration.
// 20200618  G4CMP-212: Add optional parameter for below-bandgap phonons
//		to be absorbed in the superconducting film.
// 20200626  G4CMP-215: Add function to encapsulate below-bandgap absorption,
//		apply to initial reflection condition for low-energy phonons.
// 20200626  G4CMP-216: In CalcQPEnergies, check for below-bandgap phonons,
//		and save them on phonon list for later re-emission.
// 20200627  In *EnergyRand(), move PDF expressions to functions; eliminate
//		mutation of E argument in PhononEnergyRand().
// 20200629  G4CMP-217: QPs below lowQPLimit should radiate phonon energy
//		down to gapEnergy before absorption.  Encapsulate this in
//		a function.
// 20201109  Add diagnostic text file (like downconversion and Luke).
// 20220928  G4CMP-323: Add bandgap of secondary absorber (quasiparticle trap)
// 		Drop requirement for material properties table at runtime.
// 20221006  G4CMP-330: Add temperature parameter with setter.
// 20221102  G4CMP-314: Add energy dependent efficiency for QP absorption.
// 20221116  G4CMP-343: Phonons which don't escape should be killed (dropped).
//		Subgap phonons should use MFP escape probability.

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


// Class constructor and destructor

G4CMPKaplanQP::G4CMPKaplanQP(G4MaterialPropertiesTable* prop, G4int vb)
  : verboseLevel(vb), filmProperties(0), filmThickness(0.), gapEnergy(0.),
    lowQPLimit(3.), subgapAbsorption(0.), absorberGap(0.),
    absorberEff(0.), absorberEffSlope(0.), phononLifetime(0.), 
    phononLifetimeSlope(0.), vSound(0.), temperature(0.) {
  if (prop) SetFilmProperties(prop);
}

G4CMPKaplanQP::~G4CMPKaplanQP() {
#ifdef G4CMP_DEBUG
  if (output.is_open()) output.close();
#endif
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
    phononLifetime =      prop->GetConstProperty("phononLifetime");
    phononLifetimeSlope = prop->GetConstProperty("phononLifetimeSlope");
    vSound =              prop->GetConstProperty("vSound");

    absorberEff =      (prop->ConstPropertyExists("absorberEff")
			  ? prop->GetConstProperty("absorberEff") : 1.);

    absorberEffSlope = (prop->ConstPropertyExists("absorberEffSlope")
			  ? prop->GetConstProperty("absorberEffSlope"): 0.);

    lowQPLimit =       (prop->ConstPropertyExists("lowQPLimit")
			? prop->GetConstProperty("lowQPLimit") : 3.);

    subgapAbsorption = (prop->ConstPropertyExists("subgapAbsorption")
			? prop->GetConstProperty("subgapAbsorption") : 0.);

    absorberGap =      (prop->ConstPropertyExists("absorberGap")
			? prop->GetConstProperty("absorberGap") : 0.);

    temperature =      (prop->ConstPropertyExists("temperature")
			? prop->GetConstProperty("temperature") : 0.);

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
  // FIXME: Make the properties table optional, but check if all data filled
  if (!ParamsReady()) {
    G4Exception("G4CMPKaplanQP::AbsorbPhonon()", "G4CMP001",
                RunMustBeAborted, "Thin film parameters not properly set.");
  }

  if (verboseLevel)
    G4cout << "G4CMPKaplanQP::AbsorbPhonon " << energy << G4endl;

  if (reflectedEnergies.size() > 0) {
    G4Exception("G4CMPKaplanQP::AbsorbPhonon", "G4CMP007", JustWarning,
                "Passed a nonempty reflectedEnergies vector.");
    // FIXME: Should we discard previous contents?
  }

#ifdef G4CMP_DEBUG
  if (!output.good()) {
    output.open("kaplanqp_stats");
    if (output.good()) {
      output << "Incident Energy [eV],Absorbed Energy [eV],"
	     << "Reflected Energy [eV],Reflected Phonons" << std::endl;
    }
  }
#endif

  // For the phonon to not break a Cooper pair, it must go 2*thickness,
  // assuming it goes exactly along the thickness direction, which is an
  // approximation.
  G4double frac = 2.0;

  // If phonon is not absorbed, reflect it back with no deposition
  if (IsSubgap(energy)) {
    G4double EDep = CalcSubgapAbsorption(energy, reflectedEnergies);
    if (EDep>0.) ReportAbsorption(energy, EDep, reflectedEnergies);

    return EDep;
  } else if (G4UniformRand() <= CalcEscapeProbability(energy, frac)) {
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

  ReportAbsorption(energy, EDep, reflectedEnergies);

  return EDep;
}

void G4CMPKaplanQP::
ReportAbsorption(G4double energy, G4double EDep,
		 const std::vector<G4double>& reflectedEnergies) const {
  G4double ERefl = std::accumulate(reflectedEnergies.begin(),
				   reflectedEnergies.end(), 0.);

#ifdef G4CMP_DEBUG
  if (output.good()) {
    output << energy/eV << "," << EDep/eV << "," << ERefl/eV << ","
	   << reflectedEnergies.size() << std::endl;
  }
#endif
  
  if (verboseLevel>1) {
    G4cout << " Phonon " << energy/eV << " deposited " << EDep/eV
	   << " reflected " << ERefl/eV << " as " << reflectedEnergies.size()
	   << " new phonons " << (energy-ERefl-EDep)/eV << " eV lost"
	   << G4endl;
  }
}

// Compute the probability of phonon reentering the crystal without breaking
// any Cooper pairs.

G4double G4CMPKaplanQP::CalcEscapeProbability(G4double energy,
					      G4double thicknessFrac) const {
  if (verboseLevel>1) {
    G4cout << "G4CMPKaplanQP::CalcEscapeProbability E " << energy
	   << " thickFrac " << thicknessFrac << G4endl;
  }

  if (gapEnergy <= 0.) return 1.;	// Skip phonons which can't be absorbed

  // Compute energy-dependent mean free path for phonons in film
  G4double mfp = vSound * phononLifetime /
                 (1. + phononLifetimeSlope * (energy/gapEnergy - 2.));

  if (verboseLevel>2) {
    G4cout << " mfp " << mfp << " returning "
	   << std::exp(-2.* thicknessFrac * filmThickness/mfp) << G4endl;
  }

  return std::exp(-2.* thicknessFrac * filmThickness/mfp);
}


// Model the phonons (phonEnergies) breaking Cooper pairs into quasiparticles
// (qpEnergies).

G4double 
G4CMPKaplanQP::CalcQPEnergies(std::vector<G4double>& phonEnergies,
			      std::vector<G4double>& qpEnergies) const {
  if (verboseLevel>1) {
    G4cout << "G4CMPKaplanQP::CalcQPEnergies QPcut " << lowQPLimit*gapEnergy
	   << G4endl;
  }

  // Phonons above the bandgap give all of its energy to the qp pair it breaks.
  G4double EDep = 0.;
  std::vector<G4double> newPhonEnergies;

  for (const G4double& E: phonEnergies) {
    if (IsSubgap(E)) {
      if (verboseLevel>2) G4cout << " Skipping phononE " << E << G4endl;
      newPhonEnergies.push_back(E);
      continue;
    }

    G4double qpE = QPEnergyRand(E);
    if (verboseLevel>2) G4cout << " phononE " << E << " qpE " << qpE << G4endl;

    EDep += CalcQPAbsorption(qpE, newPhonEnergies, qpEnergies);
    EDep += CalcQPAbsorption(E-qpE, newPhonEnergies, qpEnergies);
  }	// for (E: ...)

  if (verboseLevel>1)
    G4cout << " replacing phonEnergies, returning EDep " << EDep << G4endl;

  phonEnergies.swap(newPhonEnergies);
  return EDep;
}


// Model the quasiparticles (qpEnergies) emitting phonons (phonEnergies) in
// the superconductor.

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
  for (const G4double& E: qpEnergies) {
    if (verboseLevel>2) G4cout << " qpE " << E;		// Report before change

    G4double phonE = PhononEnergyRand(E);
    G4double qpE = E - phonE;
    if (verboseLevel>2)
      G4cout << " phononE " << phonE << " qpE " << qpE << G4endl;

    if (IsSubgap(phonE)) {
      EDep += CalcSubgapAbsorption(phonE, phonEnergies);
    } else {
      if (verboseLevel>2) G4cout << " Store phonE in phonEnergies" << G4endl;
      phonEnergies.push_back(phonE);
    }

    EDep += CalcQPAbsorption(qpE, phonEnergies, newQPEnergies);
  }	// for (E: ...)

  if (verboseLevel>1)
    G4cout << " replacing qpEnergies, returning EDep " << EDep << G4endl;

  qpEnergies.swap(newQPEnergies);
  return EDep;
}


// Calculate energies of phonon tracks that have reentered the crystal.

void G4CMPKaplanQP::
CalcReflectedPhononEnergies(std::vector<G4double>& phonEnergies,
                            std::vector<G4double>& reflectedEnergies) const {
  if (verboseLevel>1)
    G4cout << "G4CMPKaplanQP::CalcReflectedPhononEnergies " << G4endl;

  // There is a 50% chance that a phonon is headed away from (toward) substrate
  std::vector<G4double> newPhonEnergies;
  for (const G4double& E: phonEnergies) {
    if (verboseLevel>2) G4cout << " phononE " << E << G4endl;

    // frac = 1.5 for phonons headed away from the subst. 0.5 for toward.
    // This assumes that, on average, the phonons are spawned at the center
    // of the superconductor, which is likely not true.
    G4double frac = (G4UniformRand() < 0.5 ? 0.5 : 1.5);
    if (G4UniformRand() < CalcEscapeProbability(E, frac)) {
      if (verboseLevel>2) G4cout << " phononE got reflected" << G4endl;
      reflectedEnergies.push_back(E);
    } else if (!IsSubgap(E)) {
      newPhonEnergies.push_back(E);
    }
  }	// for (E: ...)

  phonEnergies.swap(newPhonEnergies);
}


// Compute probability of absorbing phonon below Cooper-pair breaking

G4double 
G4CMPKaplanQP::CalcSubgapAbsorption(G4double energy,
				    std::vector<G4double>& keepEnergies) const {
  if (energy < 2.*absorberGap) {	// Below absorber should just be killed
    if (verboseLevel>2)
      G4cout << " Kill phonon " << energy << " below absorber gap" << G4endl;
    
    return 0.;
  }
  
  if (G4UniformRand() < subgapAbsorption) {
    if (verboseLevel>2)
      G4cout << " Deposit phonon " << energy << " as heat" << G4endl;

    return energy;
  }
  
  if (verboseLevel>2)
    G4cout << " Record phonon " << energy << " for processing" << G4endl;

  keepEnergies.push_back(energy);
  return 0.;
}


// Handle absorption of quasiparticle energies below Cooper-pair breaking
// If qpEnergy < 3*Delta, radiate a phonon, absorb bandgap minimum

G4double 
G4CMPKaplanQP::CalcQPAbsorption(G4double qpE,
				std::vector<G4double>& phonEnergies,
				std::vector<G4double>& qpEnergies) const {
  G4double EDep = 0.;		// Energy lost by this QP into the film
  
  if (G4UniformRand() > CalcQPEfficiency(qpE)) return 0.;

  if (qpE >= lowQPLimit*gapEnergy) {
    if (verboseLevel>2) G4cout << " Storing qpE in qpEnergies" << G4endl;
    qpEnergies.push_back(qpE);
  } else if (qpE > gapEnergy) {
    if (verboseLevel>2) G4cout << " Reducing qpE to gapEnergy" << G4endl;
    EDep += CalcSubgapAbsorption(qpE-gapEnergy, phonEnergies);
    EDep += gapEnergy;
  } else {
    EDep += qpE;
  }

  return EDep;
}


// Handle quasiparticle energy-dependent absorption efficiency

G4double G4CMPKaplanQP::CalcQPEfficiency(G4double qpE) const {
    G4double eff = absorberEff + absorberEffSlope * qpE/gapEnergy;

  if (verboseLevel>2) {
    G4cout << " CalcQPEfficiency " << eff << G4endl;
  }

  return std::max(0., std::min(eff, 1.));
}


// Compute quasiparticle energy distribution from broken Cooper pair.

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
  G4double ymax = QPEnergyPDF(Energy, xmin);

  G4double xtest=0., ytest=ymax;
  do {
    ytest = G4UniformRand()*ymax;
    xtest = G4UniformRand()*(xmax-xmin) + xmin;
  } while (ytest > QPEnergyPDF(Energy, xtest));

  return xtest;
}

G4double G4CMPKaplanQP::QPEnergyPDF(G4double E, G4double x) const {
  const G4double gapsq = gapEnergy*gapEnergy;
  return ( (x*(E-x) + gapsq) / sqrt((x*x-gapsq) * ((E-x)*(E-x)-gapsq)) );
}


// Compute phonon energy distribution from quasiparticle in superconductor.

G4double G4CMPKaplanQP::PhononEnergyRand(G4double Energy) const {
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
  G4double ymax = PhononEnergyPDF(Energy, xmin);

  G4double xtest=0., ytest=ymax;
  do {
    ytest = G4UniformRand()*ymax;
    xtest = G4UniformRand()*(xmax-xmin) + xmin;
  } while (ytest > PhononEnergyPDF(Energy, xtest));

  return Energy-xtest;
}

G4double G4CMPKaplanQP::PhononEnergyPDF(G4double E, G4double x) const {
  const G4double gapsq = gapEnergy*gapEnergy;
  return ( x*(E-x)*(E-x) * (x-gapsq/E) / sqrt(x*x - gapsq) );
}

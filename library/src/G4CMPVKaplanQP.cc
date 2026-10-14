/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPVKaplanQP.cc
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
/// | absorberEff   | QP absorption efficiency     | 0.3          |
/// | absorberEffSlope    | Efficiency vs. energy        | 0.                 |
/// | temperature         | Temperature of film          | 0.05e-3*K          |
//
// $Id$
//
//
// 20260710  G4CMP-647 -- Implement virtual base class for phonon-qp
// interactions

#include "globals.hh"
#include "G4CMPVKaplanQP.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPUtils.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "Randomize.hh"
#include <numeric>

// Class constructor and destructor

G4CMPVKaplanQP::G4CMPVKaplanQP(G4MaterialPropertiesTable* prop, G4int vb)
  : verboseLevel(vb), keepAllPhonons(true),
    filmProperties(0), filmThickness(0.), gapEnergy(0.),
    lowQPLimit(3.), highQPLimit(0.), directAbsorption(0.), absorberGap(0.),
    absorberEff(1.), absorberEffSlope(0.), phononLifetime(0.), 
    phononLifetimeSlope(0.), vSound(0.), temperature(0.) {
  if (prop) SetFilmProperties(prop);
}

G4CMPVKaplanQP::~G4CMPVKaplanQP() {
#ifdef G4CMP_DEBUG
  if (output.is_open()) output.close();
#endif
}

// Configure thin film (QET, metalization, etc.) for phonon absorption

void G4CMPVKaplanQP::SetFilmProperties(G4MaterialPropertiesTable* prop) {
  if (!prop) {
    G4Exception("G4CMPVKaplanQP::SetFilmProperties()", "G4CMP001",
                RunMustBeAborted, "Null MaterialPropertiesTable vector.");
  }

  // Check that the MaterialPropertiesTable has everything we need. If it came
  // from a G4CMPSurfaceProperty, then it will be fine.
  if (!(prop->ConstPropertyExists("gapEnergy") &&
        prop->ConstPropertyExists("phononLifetime") &&
        prop->ConstPropertyExists("phononLifetimeSlope") &&
        prop->ConstPropertyExists("vSound") &&
        prop->ConstPropertyExists("filmThickness"))) {
    G4Exception("G4CMPVKaplanQP::SetFilmProperties()", "G4CMP002",
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

    highQPLimit =      (prop->ConstPropertyExists("highQPLimit")
      ? prop->GetConstProperty("highQPLimit") : 0.);

    // Backward compatible -- support both old "subgap" and new "direct" names
    directAbsorption = (prop->ConstPropertyExists("directAbsorption")
      ? prop->GetConstProperty("directAbsorption")
      : (prop->ConstPropertyExists("subgapAbsorption")
         ? prop->GetConstProperty("subgapAbsorption") : 0.)
      );

    absorberGap =      (prop->ConstPropertyExists("absorberGap")
      ? prop->GetConstProperty("absorberGap") : 0.);

    temperature =      (prop->ConstPropertyExists("temperature")
      ? prop->GetConstProperty("temperature")
      : G4CMPConfigManager::GetTemperature() );

    filmProperties = prop;
  }
}

void G4CMPVKaplanQP::
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

  G4double delta = energy-ERefl-EDep;
  if (fabs(delta) < 1e-20) delta = 0.;  // Suppress floating-point fluctuation

  if (verboseLevel>1) {
    G4cout << " Phonon " << energy/eV << " deposited " << EDep/eV
     << " reflected " << ERefl/eV << " as " << reflectedEnergies.size()
     << " new phonons " << delta/eV << " eV lost"
     << G4endl;
  }

  if (delta < 0.) {   // Actual energy excess
    G4cerr << "WARNING G4CMPVKaplanQP has excess " << delta/eV << " eV"
     << " above incident phonon." << G4endl;
  }
}

// Compute quasiparticle energy distribution from broken Cooper pair.

G4double G4CMPVKaplanQP::QPEnergyRand(G4double Energy) const {
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

G4double G4CMPVKaplanQP::QPEnergyPDF(G4double E, G4double x) const {
  const G4double gapsq = gapEnergy*gapEnergy;
  const G4double mbE = G4CMP::MaxwellBoltzmannPDF(temperature, E);
  const G4double mbEx = G4CMP::MaxwellBoltzmannPDF(temperature, E-x);
  const G4double occupy = 1. - mbE - mbEx;

  return ( occupy * (x*(E-x)+gapsq) / sqrt((x*x-gapsq) * ((E-x)*(E-x)-gapsq)) );
}


// Compute phonon energy distribution from quasiparticle in superconductor.
// NOTE:  Technically, this is the energy of the QP after emission; the
//        phonon's own energy is Ephonon = Energy - E', below

G4double G4CMPVKaplanQP::PhononEnergyRand(G4double Energy) const {
  // PDF is not integrable, so we can't do an inverse transform sampling.
  // Instead, we'll do a rejection method.
  //
  // PDF(E') = ((Energy-E')*(Energy-E') * (E'-gapEnergy*gapEnergy/Energy))
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

G4double G4CMPVKaplanQP::PhononEnergyPDF(G4double E, G4double x) const {
  const G4double gapsq = gapEnergy*gapEnergy;
  return ( (E-x)*(E-x) * (x-gapsq/E) / sqrt(x*x - gapsq) );
}

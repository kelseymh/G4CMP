/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPEnergyPartition.cc
/// \brief Implementation of the G4CMPEnergyPartition class
///   Functionality to convert energy deposition from Geant4 (both total
///   and non-ionizing) into phonons and charge carrier pairs
///
// $Id$
//
// 20160830  Fix 'A' parameter in Lindhard to convert from g/mole units.

#include "G4CMPEnergyPartition.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4DynamicParticle.hh"
#include "G4LatticePhysical.hh"
#include "G4Material.hh"
#include "G4PhononPolarization.hh"
#include "G4Pow.hh"
#include "G4PrimaryParticle.hh"
#include "G4RandomDirection.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <cmath>
#include <vector>


// Constructors and destructor

G4CMPEnergyPartition::G4CMPEnergyPartition(G4Material* mat,
					   G4LatticePhysical* lat)
  : G4CMPProcessUtils(), material(mat), holeFraction(0.5),
    verboseLevel(G4CMPConfigManager::GetVerboseLevel()),
    nPairs(0), chargeEnergyLeft(0.), nPhonons(0), phononEnergyLeft(0.) {
  SetLattice(lat);
}

G4CMPEnergyPartition::~G4CMPEnergyPartition() {;}


// Fraction of total energy deposit in material which goes to e/h pairs

G4double G4CMPEnergyPartition::LindhardScalingFactor(G4double E) const {
  if (!material) {
    static G4bool report = true;
    if (report) {
      G4cerr << "G4CMPEnergyPartition: No material configured" << G4endl;
      report = false;
    }
    return 1.;
  }

  static const G4Pow* g4pow = G4Pow::GetInstance();	// Tabulated for speed
  const G4double Z=material->GetZ(), A=material->GetA()/(g/mole);

  if (verboseLevel>1) {
    G4cout << " LindhardScalingFactor " << E << " for (Z,A) " << Z
	   << " " << A << G4endl;
  }

  // From Lewin and Smith, 1996
  G4double epsilon = 0.0115 * E * g4pow->powN(g4pow->A13(Z),-7);
  G4double k = 0.133 * g4pow->A23(Z) / std::sqrt(A);
  G4double h = (0.7*g4pow->powA(epsilon,0.6) + 3.*g4pow->powA(epsilon,0.15)
		+ epsilon);

  if (verboseLevel>2) {
    G4cout << " eps " << epsilon << " k " << k << " h " << h << " : "
	   << 1.-1./(k*h) << G4endl;
  }

  return (1. - 1./(k*h));
}


// Apply Fano factor to convert true energy deposition to random pairs

G4double G4CMPEnergyPartition::MeasuredChargeEnergy(G4double eTrue) const {
  // Fano noise changes the measured charge energy
  // Std deviation of energy distribution
  G4double sigmaE = std::sqrt(eTrue * theLattice->GetFanoFactor()
			      * theLattice->GetPairProductionEnergy());
  return G4RandGauss::shoot(eTrue, sigmaE);
}


// Generate charge carriers and phonons, depending on interaction type

void G4CMPEnergyPartition::DoPartition(G4int PDGcode, G4double energy,
				       G4double eNIEL) {
  if (verboseLevel) {
    G4cout << "G4CMPEnergyPartition::DoPartition " << PDGcode
	   << " eTotal " << energy/MeV << " eNIEL " << eNIEL/MeV << " MeV"
	   << G4endl;
  }

  // User specified phonon energy directly; assume it is correct
  if (eNIEL > 0.) DoPartition(energy-eNIEL, eNIEL);
  else {
    if (PDGcode == 2112 || PDGcode > 10000) {	// Neutron or nucleus
      if (verboseLevel>1)
	G4cout << " Nuclear Recoil: type = " << PDGcode << G4endl;
      NuclearRecoil(energy);
    } else {
      Ionization(energy);
    }
  }
}


// Generate charge carriers and phonons according to uniform phase space

void G4CMPEnergyPartition::DoPartition(G4double eIon, G4double eNIEL) {
  if (verboseLevel>1) {
    G4cout << "G4CMPEnergyPartition::DoPartition eIon " << eIon/MeV
	   << " eNIEL " << eNIEL/MeV << " MeV" << G4endl;
  }

  particles.clear();		// Discard previous results

  GenerateCharges(eIon);
  GeneratePhonons(eNIEL + chargeEnergyLeft);
}

void G4CMPEnergyPartition::GenerateCharges(G4double energy) {
  if (verboseLevel) G4cout << " GenerateCharges " << energy << G4endl;

  G4double ePair = theLattice->GetPairProductionEnergy();
  G4double eMeas = MeasuredChargeEnergy(energy);	// Applies Fano factor

  nPairs = std::floor(eMeas / ePair);		// Average number of e/h pairs
  particles.reserve(particles.size() + 2*nPairs);

  if (verboseLevel>1)
    G4cout << " eMeas " << eMeas/MeV << " MeV => " << nPairs << " pairs"
	   << G4endl;

  chargeEnergyLeft = eMeas;
  while (chargeEnergyLeft > ePair) {
    AddChargePair(ePair);
    chargeEnergyLeft -= ePair;
  }

  if (verboseLevel>1) G4cout << " " << chargeEnergyLeft << " excess" << G4endl;
}

void G4CMPEnergyPartition::AddChargePair(G4double ePair) {
  G4double eFree = ePair - theLattice->GetBandGapEnergy(); // TODO: Is this right?

  particles.push_back(Data(G4CMPDriftElectron::Definition(),G4RandomDirection(),
			 (1.-holeFraction)*eFree));

  particles.push_back(Data(G4CMPDriftHole::Definition(), G4RandomDirection(),
			 holeFraction*eFree));
}

void G4CMPEnergyPartition::GeneratePhonons(G4double energy) {
  if (verboseLevel) G4cout << " GeneratePhonons " << energy << G4endl;

  G4double ePhon = theLattice->GetDebyeEnergy(); // TODO: No fluctuations yet!

  nPhonons = std::ceil(energy / ePhon);		// Average number of phonons
  particles.reserve(particles.size() + nPhonons);

  if (verboseLevel>1)
    G4cout << " ePhon " << ePhon/eV << " eV => " << nPhonons << " phonons"
	   << G4endl;

  phononEnergyLeft = energy;
  while (phononEnergyLeft > ePhon) {
    AddPhonon(ePhon);
    phononEnergyLeft -= ePhon;
  }

  if (phononEnergyLeft > 0.) AddPhonon(phononEnergyLeft);	// Residual
}

void G4CMPEnergyPartition::AddPhonon(G4double ePhon) {
  G4ParticleDefinition* pd =
    G4PhononPolarization::Get(ChoosePhononPolarization());

  particles.push_back(Data(pd, G4RandomDirection(), ePhon));
}


// Return either primary or secondary particles from partitioning

void G4CMPEnergyPartition::
GetPrimaries(std::vector<G4PrimaryParticle*>& primaries) const {
  if (verboseLevel) G4cout << "G4CMPEnergyPartition::GetPrimaries" << G4endl;

  primaries.clear();
  primaries.resize(particles.size());	// Pre-create list, address by index

  if (verboseLevel>1) G4cout << " filling " << particles.size() << G4endl;

  for (size_t i=0; i<particles.size(); i++) {
    const Data& p = particles[i];	// For convenience below
    primaries[i] = new G4PrimaryParticle();
    primaries[i]->SetParticleDefinition(p.pd);
    primaries[i]->SetMomentumDirection(p.dir);
    primaries[i]->SetKineticEnergy(p.ekin);

    if (verboseLevel==3) {
      G4cout << i << " : " << p.pd->GetParticleName() << " " << p.ekin/eV
	     << " eV along " << p.dir << G4endl;
    } else if (verboseLevel>3) {
      G4cout << i << " : ";
      primaries[i]->Print();
    }
  }
}

void G4CMPEnergyPartition::
GetSecondaries(std::vector<G4Track*>& secondaries) const {
  if (verboseLevel) G4cout << "G4CMPEnergyPartition::GetSecondaries" << G4endl;

  secondaries.clear();
  secondaries.resize(particles.size());	// Pre-create list, address by index

  if (verboseLevel>1) G4cout << " filling " << particles.size() << G4endl;

  for (size_t i=0; i<particles.size(); i++) {
    const Data& p = particles[i];	// For convenience below
    secondaries[i] = CreateTrack(p.pd, p.dir, p.ekin);

    if (verboseLevel==3) {
      G4cout << i << " : " << p.pd->GetParticleName() << " " << p.ekin/eV
	     << " eV along " << p.dir << G4endl;
    } else if (verboseLevel>3) {
      G4cout << i << " : ";
      secondaries[i]->GetDynamicParticle()->DumpInfo();
    }
  }
}

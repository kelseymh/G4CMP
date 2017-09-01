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
// 20160830  Apply production biasing for primaries and secondaries
// 20160830  Fix 'A' parameter in Lindhard to convert from g/mole units.
// 20170524  Add constructor and accessor for position argument
// 20170728  Forgot to assign material to data member in ctor.
// 20170731  Move point-to-volume conversion to G4CMPGeometryUtils.
// 20170802  Add constructor and accessor for volume argument, particle change
// 20170830  Use downsampling energy scale parameter in DoPartition()
// 20170901  Add support for putting primaries directly into event

#include "G4CMPEnergyPartition.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4CMPGeometryUtils.hh"
#include "G4CMPSecondaryUtils.hh"
#include "G4CMPUtils.hh"
#include "G4DynamicParticle.hh"
#include "G4Event.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4Navigator.hh"
#include "G4PhononPolarization.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4RandomDirection.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4VParticleChange.hh"
#include "G4VPhysicalVolume.hh"
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

G4CMPEnergyPartition::G4CMPEnergyPartition(const G4VPhysicalVolume* volume)
  : G4CMPEnergyPartition() {
  UseVolume(volume);
}

G4CMPEnergyPartition::G4CMPEnergyPartition(const G4ThreeVector& pos)
  : G4CMPEnergyPartition() {
  UsePosition(pos);
}

G4CMPEnergyPartition::~G4CMPEnergyPartition() {;}


// Extract material and lattice information from geometry

void G4CMPEnergyPartition::UseVolume(const G4VPhysicalVolume* volume) {
  FindLattice(volume);
  SetMaterial(volume->GetLogicalVolume()->GetMaterial());
}

void G4CMPEnergyPartition::UsePosition(const G4ThreeVector& pos) {
  G4VPhysicalVolume* volume = G4CMP::GetVolumeAtPoint(pos);
  if (verboseLevel) 
    G4cout << "G4CMPEnergyPartition: " << pos << " in volume "
	   << volume->GetName() << G4endl;

  UseVolume(volume);
}


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

  const G4double Z=material->GetZ(), A=material->GetA()/(g/mole);

  if (verboseLevel>1) {
    G4cout << " LindhardScalingFactor " << E << " for (Z,A) " << Z
	   << " " << A << G4endl;
  }

  // From Lewin and Smith, 1996
  G4double epsilon = 0.0115 * E * std::pow(Z, -7./3.);
  G4double k = 0.133 * std::pow(Z, 2./3.) / std::sqrt(A);
  G4double h = (0.7*std::pow(epsilon,0.6) + 3.*std::pow(epsilon,0.15)
		+ epsilon);

  if (verboseLevel>2) {
    G4cout << " eps " << epsilon << " k " << k << " h " << h << " : "
	   << k*h / (1.+k*h) << G4endl;
  }

  return (k*h / (1.+k*h));
}


// Apply Fano factor to convert true energy deposition to random pairs

G4double G4CMPEnergyPartition::MeasuredChargeEnergy(G4double eTrue) const {
  // Fano noise changes the measured charge energy
  // Std deviation of energy distribution

  if (!G4CMPConfigManager::FanoStatisticsEnabled()) {
    return eTrue;
  }

  G4double sigmaE = std::sqrt(eTrue * theLattice->GetFanoFactor()
                              * theLattice->GetPairProductionEnergy());
  return G4RandGauss::shoot(eTrue, sigmaE);
}


// Generate charge carriers and phonons, depending on interaction type

void G4CMPEnergyPartition::DoPartition(G4int PDGcode, G4double energy,
				       G4double eNIEL) {
  if (verboseLevel) {
    G4cout << "G4CMPEnergyPartition::DoPartition: ParticleID " << PDGcode
     << "; eTotal " << energy/MeV << " MeV; eNIEL " << eNIEL/MeV << " MeV"
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
  G4double samplingScale = G4CMPConfigManager::GetSamplingEnergy();

  if (verboseLevel>1) {
    G4cout << "G4CMPEnergyPartition::DoPartition: eIon " << eIon/MeV
	   << " MeV, eNIEL " << eNIEL/MeV << " MeV" << G4endl;
  }

  particles.clear();		// Discard previous results

  // Apply downsampling if total energy is above scale
  if (samplingScale > 0.) ComputeDownsampling(eIon, eNIEL);

  GenerateCharges(eIon);
  GeneratePhonons(eNIEL + chargeEnergyLeft);
}

void G4CMPEnergyPartition::ComputeDownsampling(G4double eIon, G4double eNIEL) {
  G4double samplingScale = G4CMPConfigManager::GetSamplingEnergy();
  if (samplingScale <= 0.) return;		// Avoid unnecessary work

  if (verboseLevel>1) {
    G4cout << "G4CMPEnergyPartition::ComputeDownsampling: scale "
	   << samplingScale/eV << " eV" << G4endl;
  }

  G4double phononSamp = (eNIEL>samplingScale) ? samplingScale/eNIEL : 1.;
  if (verboseLevel>2)
    G4cout << " Downsample " << phononSamp << " for primary phonons" << G4endl;

  G4CMPConfigManager::SetGenPhonons(phononSamp);

  G4double chargeSamp = (eIon>samplingScale)? samplingScale/eIon : 1.;
  if (verboseLevel>2)
    G4cout << " Downsample " << chargeSamp << " for primary charges" << G4endl;

  G4CMPConfigManager::SetGenCharges(chargeSamp);

  // FIXME:  Want to estimate # Luke phonons per charge carrier
  G4double lukeSamp = chargeSamp;
  if (verboseLevel>2)
    G4cout << " Downsample " << lukeSamp << " Luke-phonon emission" << G4endl;

  G4CMPConfigManager::SetLukeSampling(lukeSamp);
}

void G4CMPEnergyPartition::GenerateCharges(G4double energy) {
  if (verboseLevel) G4cout << " GenerateCharges " << energy/MeV << " MeV" << G4endl;

  G4double ePair = theLattice->GetPairProductionEnergy();
  G4double eMeas = MeasuredChargeEnergy(energy);	// Applies Fano factor

  if (eMeas > 0) {
    nPairs = std::floor(eMeas / ePair);   // Average number of e/h pairs
  } else {
    nPairs = 0; // This prevents nPairs from blowing up when assigned as a negative number
                // which would cause particles.reserve() to crash
  }
  
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
  if (verboseLevel) G4cout << " GeneratePhonons " << energy/MeV << " MeV" <<  G4endl;

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


// Return primary particles from partitioning as list

void G4CMPEnergyPartition::
GetPrimaries(std::vector<G4PrimaryParticle*>& primaries) const {
  if (verboseLevel) G4cout << "G4CMPEnergyPartition::GetPrimaries" << G4endl;

  primaries.clear();
  primaries.reserve(particles.size());

  if (verboseLevel>1) {
    G4cout << " processing " << particles.size()
	   << " bias " << G4CMPConfigManager::GetGenPhonons() << " phonons"
	   << " bias " << G4CMPConfigManager::GetGenCharges() << " e/h pairs"
	   << G4endl;
  }

  G4double weight = 0.;
  G4PrimaryParticle* thePrim = 0;
  for (size_t i=0; i<particles.size(); i++) {
    const Data& p = particles[i];	// For convenience below

    weight = G4CMP::ChooseWeight(p.pd);
    if (weight == 0.) continue; // Biasing rejected particle creation

    thePrim = new G4PrimaryParticle();
    thePrim->SetParticleDefinition(p.pd);
    thePrim->SetMomentumDirection(p.dir);
    thePrim->SetKineticEnergy(p.ekin);
    thePrim->SetWeight(weight);
    primaries.push_back(thePrim);

    if (verboseLevel==3) {
      G4cout << i << " : " << p.pd->GetParticleName() << " " << p.ekin/eV
	           << " eV along " << p.dir << " (w " << weight << ")" << G4endl;
    } else if (verboseLevel>3) {
      G4cout << i << " : ";
      thePrim->Print();
    }
  }

  primaries.shrink_to_fit();		// Reduce footprint if biasing done
}

// Return primary particles from partitioning directly into event

void G4CMPEnergyPartition::
GetPrimaries(G4Event* event, const G4ThreeVector& pos, G4double time,
	     G4int maxPerVertex) const {
  if (verboseLevel) {
    G4cout << "G4CMPEnergyPartition::GetPrimaries @ " << pos
	   << " up to " << maxPerVertex << " per vertex" << G4endl;
  }

  std::vector<G4PrimaryParticle*> primaries;	// Can we make this mutable?
  GetPrimaries(primaries);

  G4double chargeEtot = 0.;		// Cumulative buffers for diagnostics
  G4double phononEtot = 0.;
  G4double bandgap = GetLattice()->GetBandGapEnergy()/2.;

  G4PrimaryVertex* vertex = CreateVertex(event, pos, time);

  for (size_t i=0; i<primaries.size(); i++) {
    // Change vertices if current one is full
    if (maxPerVertex>0 && vertex->GetNumberOfParticle() > maxPerVertex) {
      vertex = CreateVertex(event, pos, time);
    }

    vertex->SetPrimary(primaries[i]);		// Add primary to vertex

    // Accumulate results for diagnostic output
    if (verboseLevel>2) {
      const G4ParticleDefinition* pd = primaries[i]->GetParticleDefinition();
      G4double ekin = primaries[i]->GetKineticEnergy();
      G4double weight = primaries[i]->GetWeight();
      if (G4CMP::IsPhonon(pd)) phononEtot += ekin * weight;
      if (G4CMP::IsChargeCarrier(pd)) chargeEtot += (ekin + bandgap) * weight;
    }
  }

  if (verboseLevel>2) {
    G4cout << "Energy in electron-hole pairs " << chargeEtot/keV << " keV\n"
           << "Energy in phonons " << phononEtot/keV << " keV\n"
	   << "Primary particles " << primaries.size()
	   << " in " << event->GetNumberOfPrimaryVertex() << " vertices"
           << G4endl;
  }
}

G4PrimaryVertex* 
G4CMPEnergyPartition::CreateVertex(G4Event* evt, const G4ThreeVector& pos,
				   G4double time) const {
  G4PrimaryVertex* vertex = new G4PrimaryVertex(pos, time);
  evt->AddPrimaryVertex(vertex);

  return vertex;
}


// Return secondary particles from partitioning as list

void G4CMPEnergyPartition::
GetSecondaries(std::vector<G4Track*>& secondaries, G4double trkWeight) const {
  if (verboseLevel) {
    G4cout << "G4CMPEnergyPartition::GetSecondaries, parent weight "
	   << trkWeight << G4endl;
  }

  secondaries.clear();
  secondaries.reserve(particles.size());

  if (verboseLevel>1) {
    G4cout << " processing " << particles.size()
	   << " bias " << G4CMPConfigManager::GetGenPhonons() << " phonons"
	   << " bias " << G4CMPConfigManager::GetGenCharges() << " e/h pairs"
	   << G4endl;
  }

  G4double weight = 0.;
  G4Track* theSec = 0;
  for (size_t i=0; i<particles.size(); i++) {
    const Data& p = particles[i];	// For convenience below

    weight = trkWeight * G4CMP::ChooseWeight(p.pd);
    if (weight == 0.) continue;

    theSec = G4CMP::CreateSecondary(*GetCurrentTrack(), p.pd, p.dir, p.ekin);
    theSec->SetWeight(weight);
    secondaries.push_back(theSec);

    if (verboseLevel==3) {
      G4cout << i << " : " << p.pd->GetParticleName() << " " << p.ekin/eV
	           << " eV along " << p.dir << " (w " << weight << ")" << G4endl;
    } else if (verboseLevel>3) {
      G4cout << i << " : ";
      theSec->GetDynamicParticle()->DumpInfo();
    }
  }

  secondaries.shrink_to_fit();		// Reduce footprint if biasing done
}

// Return secondary particles from partitioning directly into event

void G4CMPEnergyPartition::
GetSecondaries(G4VParticleChange* aParticleChange) const {
  if (verboseLevel) {
    G4cout << "G4CMPEnergyPartition::GetSecondaries into ParticleChange"
	   << G4endl;
  }

  std::vector<G4Track*> secondaries;	// Can we make this a mutable buffer?
  GetSecondaries(secondaries, aParticleChange->GetWeight());

  aParticleChange->SetNumberOfSecondaries(secondaries.size());
  aParticleChange->SetSecondaryWeightByProcess(true);
  
  while (!secondaries.empty()) {	// Move entries from list to tracking
    aParticleChange->AddSecondary(secondaries.back());
    secondaries.pop_back();
  }
}


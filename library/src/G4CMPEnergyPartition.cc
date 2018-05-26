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
// 20170925  Add support for distributing charges around position
// 20171213  Apply downsampling up front, to initial generated Data objects
// 20180424  Count actual number of tracks produced during downsampling,
//		set their weights to the ratio of true/produced.
// 20180503  Protect against negative "energy left".
// 20180511  Protect GetSecondaries()/GetPrimaries() from zero generated.

#include "G4CMPEnergyPartition.hh"
#include "G4CMPChargeCloud.hh"
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
  : G4CMPProcessUtils(), verboseLevel(G4CMPConfigManager::GetVerboseLevel()),
    material(mat), holeFraction(0.5), nParticlesMinimum(10),
    cloud(new G4CMPChargeCloud), nCharges(0), nPairs(0),
    chargeEnergyLeft(0.), nPhonons(0), phononEnergyLeft(0.) {
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

G4CMPEnergyPartition::~G4CMPEnergyPartition() {
  delete cloud; cloud=0;
}


// Extract material and lattice information from geometry

void G4CMPEnergyPartition::UseVolume(const G4VPhysicalVolume* volume) {
  FindLattice(volume);
  SetMaterial(volume->GetLogicalVolume()->GetMaterial());
  cloud->UseVolume(volume);
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
  nPairs = nPhonons = 0;

  // Apply downsampling if total energy is above scale
  if (samplingScale > 0.) ComputeDownsampling(eIon, eNIEL);

  chargeEnergyLeft = 0.;
  GenerateCharges(eIon);
  GeneratePhonons(eNIEL + chargeEnergyLeft);

  particles.shrink_to_fit();	// Reduce size to match generated particles
}

void G4CMPEnergyPartition::ComputeDownsampling(G4double eIon, G4double eNIEL) {
  G4double samplingScale = G4CMPConfigManager::GetSamplingEnergy();
  if (samplingScale <= 0.) return;		// Avoid unnecessary work

  if (verboseLevel>1) {
    G4cout << "G4CMPEnergyPartition::ComputeDownsampling: scale energy to "
	   << samplingScale/eV << " eV" << G4endl;
  }

  // Compute phonon scaling factor only if not fully suppressed
  if (G4CMPConfigManager::GetGenPhonons() > 0.) {
    G4double phononSamp = (eNIEL>samplingScale) ? samplingScale/eNIEL : 1.;
    if (verboseLevel>2)
      G4cout << " Downsample " << phononSamp << " primary phonons" << G4endl;

    G4CMPConfigManager::SetGenPhonons(phononSamp);
  }

  // Compute charge scaling factor only if not fully suppressed
  if (G4CMPConfigManager::GetGenCharges() > 0.) {
    G4double chargeSamp = (eIon>samplingScale)? samplingScale/eIon : 1.;
    if (verboseLevel>2)
      G4cout << " Downsample " << chargeSamp << " primary charges" << G4endl;
    
    G4CMPConfigManager::SetGenCharges(chargeSamp);

    // FIXME:  Want to estimate # Luke phonons per charge carrier
    G4double lukeSamp = chargeSamp;
    if (verboseLevel>2)
      G4cout << " Downsample " << lukeSamp << " Luke-phonon emission" << G4endl;
    
    G4CMPConfigManager::SetLukeSampling(lukeSamp);
  }
}

void G4CMPEnergyPartition::GenerateCharges(G4double energy) {
  if (G4CMPConfigManager::GetGenCharges() <= 0.) return;	// Suppressed

  if (verboseLevel)
    G4cout << " GenerateCharges " << energy/MeV << " MeV" << G4endl;

  G4double ePair = theLattice->GetPairProductionEnergy();
  G4double eMeas = MeasuredChargeEnergy(energy);	// Applies Fano factor

  if (eMeas > 0) {
    nPairs = std::floor(eMeas / ePair);   // Average number of e/h pairs
  } else {
    eMeas = 0.;
    nPairs = 0; // This prevents nPairs from blowing up when negative
  }
  
  // Only apply downsampling to sufficiently large statistics
  G4double scale = (nPairs<=nParticlesMinimum ? 1.
		    : G4CMPConfigManager::GetGenCharges());

  if (verboseLevel>1) {
    G4cout << " eMeas " << eMeas/MeV << " MeV => " << nPairs << " pairs"
	   << " downsample " << scale << G4endl;
  }

  chargeEnergyLeft = eMeas;
  nCharges = 0;
  if (nPairs == 0) return;		// No charges could be produced

  particles.reserve(particles.size() + int(2.2*scale*nPairs));	// 10% overhead

  // For downsampling, ensure that there are sufficient charge pairs
  while (nCharges < scale*nPairs) {
    if (nCharges > 0)
      particles.erase(particles.end()-2*nCharges, particles.end());

    chargeEnergyLeft = eMeas;
    nCharges = 0;

    while (chargeEnergyLeft >= ePair) {
      if (G4UniformRand()<scale) {	// Apply downsampling up front
	AddChargePair(ePair);
	nCharges++;
      }

      chargeEnergyLeft -= ePair;
    }	// while (chargeEnergyLeft

    if (verboseLevel>2)
      G4cout << " generated " << nCharges << " e-h pairs" << G4endl;
  }	// while (nCharges==0

  if (chargeEnergyLeft < 0.) chargeEnergyLeft = 0.;	// Avoid round-offs

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
  if (G4CMPConfigManager::GetGenPhonons() <= 0.) return;	// Suppressed

  if (verboseLevel)
    G4cout << " GeneratePhonons " << energy/MeV << " MeV" <<  G4endl;

  G4double ePhon = theLattice->GetDebyeEnergy(); // TODO: No fluctuations yet!

  nPhonons = std::ceil(energy / ePhon);		// Average number of phonons
  ePhon = energy / nPhonons;			// Split energy evenly to all

  // Only apply downsampling to sufficiently large statistics
  G4double scale = (nPhonons<=nParticlesMinimum ? 1.
		    : G4CMPConfigManager::GetGenPhonons());

  if (verboseLevel>1) {
    G4cout << " ePhon " << ePhon/eV << " eV => " << nPhonons << " phonons"
	   << " downsample " << scale << G4endl;
  }

  particles.reserve(particles.size() + int(1.1*scale*nPhonons)); // 10% overhead

  // For downsampling, ensure that there are sufficient phonons
  size_t nGenPhonons = 0;
  while (nGenPhonons < scale*nPhonons/2) {
    if (nGenPhonons > 0)
      particles.erase(particles.end()-nGenPhonons, particles.end());

    phononEnergyLeft = energy;
    nGenPhonons = 0;
    
    while (phononEnergyLeft >= ePhon) {
      if (G4UniformRand()<scale) {	// Apply downsampling up front
	AddPhonon(ePhon);
	nGenPhonons++;
      }
      
      phononEnergyLeft -= ePhon;

      if (phononEnergyLeft > 0 && phononEnergyLeft < ePhon)
	phononEnergyLeft = ePhon;	// Correct for round-off issues
    }	// while (phononEnergyLeft
    
    if (verboseLevel>2)
      G4cout << " generated " << nGenPhonons << " phonons" << G4endl;
  }	// while (nGenPhonons
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

  if (verboseLevel>1) G4cout << " processing " << particles.size() << G4endl;

  // Get number of generated phonons to compute weight below
  size_t nGenPhonons = particles.size() - 2*nCharges;
  G4double phononWt = nGenPhonons>0 ? G4double(nPhonons)/nGenPhonons : 0.;
  G4double chargeWt = nCharges>0 ? G4double(nPairs)/nCharges : 0.;

  G4double weight = 0.;
  G4PrimaryParticle* thePrim = 0;
  for (size_t i=0; i<particles.size(); i++) {
    const Data& p = particles[i];	// For convenience below

    // Set weight so that generated particles map back to expected true number
    weight = (G4CMP::IsPhonon(p.pd) ? phononWt : chargeWt);

    thePrim = new G4PrimaryParticle();
    thePrim->SetParticleDefinition(p.pd);
    thePrim->SetMomentumDirection(p.dir);
    thePrim->SetKineticEnergy(p.ekin);
    thePrim->SetWeight(weight);
    primaries.push_back(thePrim);

    if (verboseLevel==3) {
      G4cout << i << " : " << p.pd->GetParticleName() << " " << p.ekin/eV
	     << " eV along " << p.dir << " (w " << weight << ")"
	     << G4endl;
    } else if (verboseLevel>3) {
      G4cout << i << " : ";
      thePrim->Print();
    }
  }
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

  // Generate charge carriers in region around track position
  G4bool doCloud = G4CMPConfigManager::CreateChargeCloud();	// Convenience
  if (doCloud) {
    cloud->SetVerboseLevel(verboseLevel);
    cloud->SetTouchable(G4CMP::CreateTouchableAtPoint(pos));
    cloud->Generate(nCharges, pos);
  }

  // Buffer for active vertices, for use with charge cloud
  std::map<G4int, G4PrimaryVertex*> activeVtx;

  G4int ichg = 0;		// Counter to track charge cloud entries
  for (size_t i=0; i<primaries.size(); i++) {
    G4bool qcloud = doCloud && !G4CMP::IsPhonon(primaries[i]->GetG4code());
    G4int chgbin = qcloud ? cloud->GetPositionBin(ichg++) : -1;

    G4PrimaryVertex*& vertex = activeVtx[chgbin];	// Ref for convenience

    // Create new vertex at pos if needed, or if current one is full
    if (!vertex ||
	(maxPerVertex>0 && vertex->GetNumberOfParticle()>maxPerVertex)) {
      G4ThreeVector binpos = chgbin>=0 ? cloud->GetBinCenter(chgbin) : pos;
      vertex = CreateVertex(event, binpos, time);
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

  // Generate charge carriers in region around track position
  G4bool doCloud = G4CMPConfigManager::CreateChargeCloud();	// Convenience
  if (doCloud) {
    cloud->SetVerboseLevel(verboseLevel);
    cloud->SetTouchable(GetCurrentTouchable());
    cloud->Generate(nCharges, GetCurrentTrack()->GetPosition());
  }

  if (verboseLevel>1) G4cout << " processing " << particles.size() << G4endl;

  // Get number of generated phonons to compute weight below
  size_t nGenPhonons = particles.size() - 2*nCharges;
  G4double phononWt = nGenPhonons>0 ? G4double(nPhonons)/nGenPhonons : 0.;
  G4double chargeWt = nCharges>0 ? G4double(nPairs)/nCharges : 0.;

  G4double weight = 0.;
  G4Track* theSec = 0;
  G4int ichg = 0;			// Index to deal with charge cloud

  for (size_t i=0; i<particles.size(); i++) {
    const Data& p = particles[i];	// For convenience below

    // Set weights so that generated particles map back to expected true number
    weight = (G4CMP::IsPhonon(p.pd) ? phononWt : chargeWt);

    theSec = G4CMP::CreateSecondary(*GetCurrentTrack(), p.pd, p.dir, p.ekin);
    theSec->SetWeight(weight);
    secondaries.push_back(theSec);

    // Adjust positions of charges according to generated distribution
    if (doCloud && G4CMP::IsChargeCarrier(theSec))
      theSec->SetPosition(cloud->GetPosition(ichg++));

    if (verboseLevel==3) {
      G4cout << i << " : " << p.pd->GetParticleName() << " " << p.ekin/eV
	     << " eV along " << p.dir << " (w " << weight << ")"
	     << G4endl;
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


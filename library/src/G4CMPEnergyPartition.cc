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
// 20180801  Add weighting bounds for computing Luke-phonon sampling.
// 20180827  Add flag to suppress use of downsampling energy scale
// 20180828  BUG FIX:  GetSecondaries() was not using trkWeight
// 20180831  Fix compiler warnings when comparing nParticlesMinimum
// 20190711  Use selectable NIEL partition function, via ConfigManager.
// 20190714  Convert PDGcode to Z and A (in amu) for use with NIEL function.
// 20191009  Produce charge pairs below pair-energy, down to bandgap.
// 20191017  Fix PDGcode usage for nuclei to look up in G4IonTable.
// 20191106  Protect against exactly zero energy passed to GeneratePhonons()
// 20200217  Fill new 'hit' container with generated parameters
// 20200219  Replace use of G4HitsCollection with singleton data container
// 20200222  Add control flag to turn off creating summary data
// 20200316  Improve calculations of charge and phonon energy summaries
// 20200328  Protect against invalid energy inputs
// 20200805  Use electric field in volume to estimate Luke gain, sampling
// 20201013  Implement Fano fluctuations as applying to Npair, not Emeas
// 20201020  Use "interpolation" to match input Fano factor and mean Npair
// 20201205  Use ApplySurfaceClearance() for primary production, to avoid
//		creating particles which escape from volume.
// 20210202  in DoPartition(PDGcode, ...) store particle type in summary.
// 20210328  Split ComputeDownsampling() into individual computation functions  
// 20210706  Add flag to control whether ComputeLukeSampling() is used.
// 20210820  Apply downsampling deterministically, not in loop over tracks,
//		rename particle count data member for clarity.  Store weight
//		for each particle in internal "Data" buffer.  Store both true
//		and downsampled particle counts in summary buffer.
// 20210820  Add estimate of NTL (Luke) emission energy to summary buffer.
// 20210915  Change ComputeLukeDownsampling() to use new parameter requesting
//		specific number of NTL phonons; estimated using nPairs.
// 20211030  Add track and step summary information to support data analysis
// 20220216  Add interface to do partitioning directly from StepAccumulator.
// 20220228  In GetSecondaries(), don't overwrite previously set position info.
// 20220818  G4CMP-309 -- Don't skip GenerateCharges() or GeneratePhonons() if
//		zero downsampling; want to get summary data filled every time.
// 20221025  G4CMP-335 -- Skip and rethrow nPairs=0 returned from FanoBinomial.

#include "G4CMPEnergyPartition.hh"
#include "G4CMPChargeCloud.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4CMPFanoBinomial.hh"
#include "G4CMPFieldUtils.hh"
#include "G4CMPGeometryUtils.hh"
#include "G4CMPPartitionData.hh"
#include "G4CMPPartitionSummary.hh"
#include "G4CMPSecondaryUtils.hh"
#include "G4CMPStepAccumulator.hh"
#include "G4CMPUtils.hh"
#include "G4VNIELPartition.hh"
#include "G4DynamicParticle.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4IonTable.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4Neutron.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PhononPolarization.hh"
#include "G4PhysicalConstants.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4RandomDirection.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VParticleChange.hh"
#include "G4VPhysicalVolume.hh"
#include "Randomize.hh"
#include "CLHEP/Random/RandBinomial.h"
#include <cmath>
#include <vector>


// Constructors and destructor

// TEMPORARY:  Flag to either compute Luke sampling, or use preset value
namespace {
  G4double lukeDownsampling = true;
}

G4CMPEnergyPartition::G4CMPEnergyPartition(G4Material* mat,
					   G4LatticePhysical* lat)
  : G4CMPProcessUtils(), verboseLevel(G4CMPConfigManager::GetVerboseLevel()),
    fillSummaryData(false), material(mat), biasVoltage(0.), 
    holeFraction(0.5), nParticlesMinimum(10),
    applyDownsampling(true), cloud(new G4CMPChargeCloud),
    nPairsTrue(0), nPairsGen(0), chargeEnergyLeft(0.),
    nPhononsTrue(0), nPhononsGen(0), phononEnergyLeft(0.),
    summary(0) {
  SetLattice(lat);

  // TEMPORARY: If user set Luke sampling negative, we compute it below
  lukeDownsampling = (G4CMPConfigManager::GetLukeSampling() < 0.);
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
  // NOTE: "summary" is owned by G4Event/G4HitsCollection, do not delete
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
    G4cout << "G4CMPEnergyPartition: " << pos << " at volume "
	   << volume->GetName() << G4endl;

  UseVolume(volume);
  SetBiasVoltage(pos);
}

void G4CMPEnergyPartition::SetBiasVoltage(const G4ThreeVector& pos) {
  G4VTouchable* touch = G4CMP::CreateTouchableAtPoint(pos);
  biasVoltage = G4CMP::GetBiasThroughPosition(touch, pos);

  if (verboseLevel) {
    G4cout << "G4CMPEnergyPartition: est. " << biasVoltage/volt << " V"
	   << " across " << touch->GetVolume()->GetName() << " @ " << pos
	   << G4endl;
  }
}


// Create summary data container and put into event "hits" collection

G4CMPPartitionData* G4CMPEnergyPartition::CreateSummary() {
  if (verboseLevel) G4cout << "G4CMPEnergyPartition::CreateSummary" << G4endl;

  if (!fillSummaryData) {		// No collection, just keep local
    if (!summary) summary = new G4CMPPartitionData;
    return summary;
  }

  summary = new G4CMPPartitionData;	// Ownership transfers to container
  G4CMPPartitionSummary::Insert(summary);

  if (verboseLevel>2) {
    G4cout << " Partition summary contains "
	   << G4CMPPartitionSummary::Entries() << " records" << G4endl;
  }

  return summary;
}


// Fraction of total energy deposit in material which goes to e/h pairs

G4double G4CMPEnergyPartition::
LindhardScalingFactor(G4double E, G4double Z, G4double A) const {
  if (!material) {
    G4Exception("G4CMPEnergyPartition", "G4CMP1000", RunMustBeAborted,
		"No material configured for energy partition");
    return 1.;
  }

  const G4VNIELPartition* nielFunc = G4CMPConfigManager::GetNIELPartition();
  return nielFunc->PartitionNIEL(E, material, Z, A);
}


// Apply Fano factor to convert true energy deposition to random pairs

G4double G4CMPEnergyPartition::MeasuredChargePairs(G4double eTrue) const {
  if (eTrue < theLattice->GetBandGapEnergy()) return 0.;
  if (eTrue <= theLattice->GetPairProductionEnergy()) return 1.;

  G4double Ntrue = eTrue/theLattice->GetPairProductionEnergy();

  // Fano noise changes the number of generated charges
  if (!G4CMPConfigManager::FanoStatisticsEnabled()) {
    summary->FanoFactor = 0.;
    return std::round(Ntrue);
  }

  // Store Fano factor from material for reference
  summary->FanoFactor = theLattice->GetFanoFactor();

  if (verboseLevel>1) {
    G4cout << "Using Ntrue " << Ntrue << " and F "  << summary->FanoFactor
	   << " for Fano binomial" << G4endl;
  }

  // Interpolated binominals to reproduce preset Fano factor
  // See https://www.slac.stanford.edu/exp/cdms/ScienceResults/DataReleases/20190401_HVeV_Run1/HVeV_R1_Data_Release_20190401.pdf
  G4double tryPairs;
  G4int maxThrow=1000;
  do {
    tryPairs = G4CMP::FanoBinomial::shoot(Ntrue, summary->FanoFactor);
  } while (tryPairs<1. && --maxThrow);
  if (tryPairs==0) tryPairs = Ntrue;		// Didn't get good throw

  return tryPairs;
}


// Generate charge carriers and phonons, depending on interaction type

void 
G4CMPEnergyPartition::DoPartition(const G4CMPStepAccumulator* steps) {
  if (!steps) return;		// Avoid unnecessary work

  if (verboseLevel) {
    G4cout << "G4CMPEnergyPartition::DoPartition steps" << G4endl
	   << *steps << G4endl;
  }

  DoPartition(steps->pd->GetPDGEncoding(), steps->Edep, steps->Eniel);

  // Store position information in summary block
  if (summary) {
    summary->position[0] = steps->end[0];
    summary->position[1] = steps->end[1];
    summary->position[2] = steps->end[2];
    summary->position[3] = steps->time;
    
    summary->trackID = steps->trackID;
    summary->stepID  = steps->stepID;
  }
}


// Generate charge carriers and phonons, depending on interaction type

void G4CMPEnergyPartition::DoPartition(G4int PDGcode, G4double energy,
				       G4double eNIEL) {
  if (verboseLevel) {
    G4cout << "G4CMPEnergyPartition::DoPartition: ParticleID " << PDGcode
     << "; eTotal " << energy/MeV << " MeV; eNIEL " << eNIEL/MeV << " MeV"
	   << G4endl;
  }

  if (energy <= 0. || eNIEL < 0.) {
    G4ExceptionDescription msg;
    msg << "Invalid energy input:";
    if (energy <= 0.) msg << " energy " << energy/eV << " eV";
    if (eNIEL < 0.) msg << " NIEL " << eNIEL/eV << " eV";

    G4Exception("G4CMPEnergyPartition::DoPartition", "Partition001",
		EventMustBeAborted, msg);
    return;
  }

  // User specified phonon energy directly; assume it is correct
  if (eNIEL > 0.) DoPartition(energy-eNIEL, eNIEL);
  else {
    G4ParticleDefinition* proj = 0;	// To get nuclear recoil info

    if (PDGcode == 2112)		// Neutron; treat as nuclear particle
      proj = G4Neutron::Definition();
    else if (PDGcode > 1000000000)	// Geant4 native nucleus encoding
      proj = G4IonTable::GetIonTable()->GetIon(PDGcode);
    else if (PDGcode > 10000)		// Nucleus pseudo-code, AAAZZZ
      proj = G4IonTable::GetIonTable()->GetIon(PDGcode%1000, PDGcode/1000);

    if (proj) {
      G4double Z=proj->GetAtomicNumber(), A=proj->GetPDGMass()/amu_c2;
      
      if (verboseLevel>1) {
        G4cout << " Nuclear Recoil: type " << PDGcode << " Z " << Z
	       << " A " << A << G4endl;
      }
      
      NuclearRecoil(energy, Z, A);
    } else {
      Ionization(energy);
    }
  }

  // After summary block created above, record particle type
  if (summary) summary->PDGcode = PDGcode;
}


// Generate charge carriers and phonons according to uniform phase space

void G4CMPEnergyPartition::DoPartition(G4double eIon, G4double eNIEL) {
  if (verboseLevel>1) {
    G4cout << "G4CMPEnergyPartition::DoPartition: eIon " << eIon/MeV
	   << " MeV, eNIEL " << eNIEL/MeV << " MeV" << G4endl;
  }

  if (eIon+eNIEL <= 0. || eIon< 0. || eNIEL < 0.) {
    G4ExceptionDescription msg;
    msg << "Invalid energy input:";
    if (eIon+eNIEL <= 0.) msg << " eTotal " << (eIon+eNIEL)/eV << " eV";
    if (eIon < 0.) msg << " eIon " << eIon/eV << " eV";
    if (eNIEL < 0.) msg << " NIEL " << eNIEL/eV << " eV";

    G4Exception("G4CMPEnergyPartition::DoPartition", "Partition001",
		EventMustBeAborted, msg);
    return;
  }

  particles.clear();		// Discard previous results
  nPairsTrue = nPhononsTrue = 0;

  // Set up summary information block in event
  CreateSummary();
  summary->totalEnergy = eIon+eNIEL;
  summary->truedEdx = eIon;
  summary->trueNIEL = eNIEL;
  summary->lindhardYield = eIon / (eIon+eNIEL);

  // Apply downsampling if requested
  if (applyDownsampling) ComputeDownsampling(eIon, eNIEL);

  summary->samplingEnergy  = G4CMPConfigManager::GetSamplingEnergy();
  summary->samplingCharges = G4CMPConfigManager::GetGenCharges();
  summary->samplingPhonons = G4CMPConfigManager::GetGenPhonons();
  summary->samplingLuke    = G4CMPConfigManager::GetLukeSampling();

  chargeEnergyLeft = eIon;
  GenerateCharges(eIon);
  GeneratePhonons(eNIEL + chargeEnergyLeft);

  particles.shrink_to_fit();	// Reduce size to match generated particles

  // Shuffle particles so they can be distributed along trajectories
  std::random_shuffle(particles.begin(), particles.end(), G4CMP::RandomIndex);

  if (verboseLevel && summary) summary->Print();
}


// Generate primary phonons from energy deposit using Lindhard scaling

void G4CMPEnergyPartition::
NuclearRecoil(G4double energy, G4double Z, G4double A) {
  G4double lind = LindhardScalingFactor(energy, Z, A);
  if (verboseLevel>1) G4cout << " Lindard partition factor " << lind << G4endl;

  DoPartition(energy*lind, energy*(1.-lind));
}


// Generate charge carriers and phonons with maximum energy scaling

void G4CMPEnergyPartition::ComputeDownsampling(G4double eIon, G4double eNIEL) {
  G4double samplingScale = G4CMPConfigManager::GetSamplingEnergy();
  if (samplingScale > 0.) {
    if (verboseLevel>1) {
      G4cout << "G4CMPEnergyPartition::ComputeDownsampling: scale energy to "
	     << samplingScale/eV << " eV" << G4endl;
    }
    
    ComputeChargeSampling(eIon);
    ComputePhononSampling(eNIEL);
    ComputeLukeSampling(eIon);
    return;
  }

  G4double maxLukeCount = G4CMPConfigManager::GetMaxLukePhonons();
  if (maxLukeCount > 0.) {
    if (verboseLevel>1) {
      G4cout << "G4CMPEnergyPartition::ComputeDownsampling: restrict Luke"
	     << " phonons to ~" << maxLukeCount << " per event" << G4endl;
    }
    
    G4double voltage = fabs(biasVoltage);
    G4double eLuke = eIon*eplus*voltage/theLattice->GetPairProductionEnergy();

    // Expect about 500 Luke phonons, ~ 2 meV each, per e/h pair per volt
    // Note: number varies with material, this estmate is best for germanium
    G4double nluke = eLuke/eV * 500;
    G4double lukeSamp = std::min(maxLukeCount/nluke, 1.);
  
    if (verboseLevel>2) {
      G4cout << " bias " << voltage << " V"
	     << " maxCount " << maxLukeCount << " phonons desired"
	     << "\n Downsample " << lukeSamp << " Luke-phonon emission"
	     << G4endl;
    }

    G4CMPConfigManager::SetLukeSampling(lukeSamp);
  }
}

// Compute phonon scaling factor only if not fully suppressed
// NOTE: Phonon sampling done to get same number as charge pairs

void
G4CMPEnergyPartition::ComputePhononSampling(G4double eNIEL) {
  G4double samplingScale = G4CMPConfigManager::GetSamplingEnergy();
  if (samplingScale <= 0.) return;		// No downsampling computation
  if (G4CMPConfigManager::GetGenPhonons() <= 0.) return;
  
  G4double phononScale = (samplingScale * theLattice->GetDebyeEnergy()
			  / theLattice->GetPairProductionEnergy());
  G4double phononSamp = (eNIEL>phononScale) ? phononScale/eNIEL : 1.;
  if (verboseLevel>2)
    G4cout << " Downsample " << phononSamp << " primary phonons" << G4endl;
  
  G4CMPConfigManager::SetGenPhonons(phononSamp);
}

// Compute charge scaling factor only if not fully suppressed

void
G4CMPEnergyPartition::ComputeChargeSampling(G4double eIon) {
  G4double samplingScale = G4CMPConfigManager::GetSamplingEnergy();
  if (samplingScale <= 0.) return;		// No downsampling computation
  if (G4CMPConfigManager::GetGenCharges() <= 0.) return;
  
  G4double chargeSamp = (eIon>samplingScale)? samplingScale/eIon : 1.;
  if (verboseLevel>2)
    G4cout << " Downsample " << chargeSamp << " primary charges" << G4endl;
  
  G4CMPConfigManager::SetGenCharges(chargeSamp);
}

// Compute Luke scaling factor only if not fully suppressed

void G4CMPEnergyPartition::ComputeLukeSampling(G4double eIon) {
  G4double samplingScale = G4CMPConfigManager::GetSamplingEnergy();
  if (samplingScale <= 0.) return;		// No downsampling computation
  if (!lukeDownsampling) return;		// User preset a fixed fraction

  // Expect about 500 Luke phonons, ~ 2 meV each, per e/h pair per volt
  // Note: number varies with material, this estmate is best for germanium
  G4double voltage = fabs(biasVoltage)/volt;
  G4double npair = ( std::min(samplingScale, eIon)
		     / theLattice->GetPairProductionEnergy() );
  G4double nluke = npair * (voltage+1.) * 500;	// <E> ~ 2 meV 

  // Scales to user-desired "maximum" (approximate) number of Luke phonons
  G4int maxCount = G4CMPConfigManager::GetMaxLukePhonons();
  if (maxCount <= 0.) maxCount = 10000.;
  G4double lukeSamp = std::min(maxCount/nluke, 1.);
  
  if (verboseLevel>2) {
    G4cout << " bias " << voltage << " V, scale " << samplingScale/eV << " eV"
	   << " maxCount " << maxCount << " phonons desired"
	   << "\n Downsample " << lukeSamp << " Luke-phonon emission" << G4endl;
  }

  G4CMPConfigManager::SetLukeSampling(lukeSamp);
}


// Divide ionization energy into electron/hole pairs, with Fano fluctuations

void G4CMPEnergyPartition::GenerateCharges(G4double energy) {
  if (verboseLevel)
    G4cout << " GenerateCharges " << energy/MeV << " MeV" << G4endl;

  G4double eBand = 1.01*theLattice->GetBandGapEnergy(); // Force visible energy
  G4double ePair = theLattice->GetPairProductionEnergy();

  // Use Fano factor to determine generated number of charge pairs
  if (energy > eBand) {
    nPairsTrue = MeasuredChargePairs(energy);	// Apply fluctuations
    ePair = energy/nPairsTrue;			// Split energy evenly to all
  } else {
    nPairsTrue = 0;
  }

  // Only apply downsampling to sufficiently large statistics
  G4double scale = G4CMPConfigManager::GetGenCharges();
  if (scale>0. && (G4int)nPairsTrue <= nParticlesMinimum) scale = 1.;

  if (verboseLevel>1) {
    G4cout << " nPairs " << nPairsTrue << " ==> ePair " << ePair/eV << " eV"
	   << " downsample " << scale << G4endl;
  }

  // Compute number of pairs to generate, adjust sampling scale to match
  nPairsGen = std::round(scale*nPairsTrue);
  scale = nPairsTrue>0 ? double(nPairsGen)/nPairsTrue : 1.;

  G4double nPairsWeighted = nPairsGen>0 ? nPairsGen/scale : 0.;

  // Create requested number of charge pairs with scaling factor
  if (nPairsGen > 0) {
    particles.reserve(particles.size() + nPairsGen);

    // Generate number of requested charge pairs, each with same energy
    for (size_t i=0; i<nPairsGen; i++) AddChargePair(ePair, 1./scale);
    
    if (verboseLevel>2)
      G4cout << " generated " << nPairsGen << " e-h pairs" << G4endl;
    
    chargeEnergyLeft = energy - ePair*nPairsWeighted;
    if (chargeEnergyLeft < 0.) chargeEnergyLeft = 0.;	// Avoid round-offs
  } else {
    chargeEnergyLeft = 0.;
  }

  if (verboseLevel>1) G4cout << " " << chargeEnergyLeft << " excess" << G4endl;

  // Store generated information in summary block
  if (summary) {
    summary->chargeEnergy = energy;
    summary->chargeFano = nPairsTrue*theLattice->GetPairProductionEnergy();
    summary->chargeGenerated = ePair*nPairsWeighted;
    summary->truePairs = nPairsTrue;
    summary->numberOfPairs = nPairsGen;
    summary->samplingCharges = scale;		// Store actual sampling used
    summary->lukeEnergyEst = nPairsWeighted * abs(biasVoltage);
  }

  // Estimate NTL (Luke) phonon emission from charge pairs
  // Assumes symmetry: each charge pair covers the full voltage bias
  if (verboseLevel>1) {
    G4cout << " estimating Luke emission for " << nPairsWeighted
	   << " e-h pairs across " << biasVoltage/volt << " V" << G4endl;
  }
}

void G4CMPEnergyPartition::AddChargePair(G4double ePair, G4double wt) {
  G4double eFree = ePair - theLattice->GetBandGapEnergy(); // TODO: Is this right?

  particles.push_back(Data(G4CMPDriftElectron::Definition(),G4RandomDirection(),
			   (1.-holeFraction)*eFree, wt));

  particles.push_back(Data(G4CMPDriftHole::Definition(), G4RandomDirection(),
			   holeFraction*eFree, wt));
}

void G4CMPEnergyPartition::GeneratePhonons(G4double energy) {
  if (energy <= 0.) {				// Avoid unnecessary work
    nPhononsTrue = nPhononsGen = 0;
    return;
  }

  if (verboseLevel)
    G4cout << " GeneratePhonons " << energy/MeV << " MeV" <<  G4endl;

  G4double ePhon = theLattice->GetDebyeEnergy(); // TODO: No fluctuations yet!

  nPhononsTrue = std::ceil(energy / ePhon);	// Average number of phonons
  ePhon = energy / nPhononsTrue;		// Split energy evenly to all

  // Only apply downsampling to sufficiently large statistics
  G4double scale = G4CMPConfigManager::GetGenPhonons();
  if (scale>0. && (G4int)nPhononsTrue <= nParticlesMinimum) scale = 1.;

  if (verboseLevel>1) {
    G4cout << " ePhon " << ePhon/eV << " eV => " << nPhononsTrue << " phonons"
	   << " downsample " << scale << G4endl;
  }

  // Compute number of phonons to generate, adjust sampling scale to match
  nPhononsGen = std::round(scale*nPhononsTrue);
  scale = nPhononsTrue>0 ? double(nPhononsGen)/nPhononsTrue : 1.;

  // Create requested number of phonons with scaling factor
  if (nPhononsGen > 0) {
    particles.reserve(particles.size() + nPhononsGen);

    // Generate number of requested charge pairs, each with same energy
    for (size_t i=0; i<nPhononsGen; i++) AddPhonon(ePhon, 1./scale);
    
    if (verboseLevel>2)
      G4cout << " generated " << nPhononsGen << " phonons" << G4endl;
  }

  // Store generated information in summary block
  if (summary) {
    summary->phononEnergy = energy;
    summary->phononGenerated = ePhon*nPhononsGen/scale;
    summary->truePhonons = nPhononsTrue;
    summary->numberOfPhonons = nPhononsGen;
    summary->samplingPhonons = scale;		// Store actual sampling used
  }
}

void G4CMPEnergyPartition::AddPhonon(G4double ePhon, G4double wt) {
  G4ParticleDefinition* pd =
    G4PhononPolarization::Get(ChoosePhononPolarization());

  particles.push_back(Data(pd, G4RandomDirection(), ePhon, wt));
}


// Return primary particles from partitioning as list

void G4CMPEnergyPartition::
GetPrimaries(std::vector<G4PrimaryParticle*>& primaries) const {
  if (verboseLevel) G4cout << "G4CMPEnergyPartition::GetPrimaries" << G4endl;

  primaries.clear();
  primaries.reserve(particles.size());

  if (verboseLevel>1) G4cout << " processing " << particles.size() << G4endl;

  G4PrimaryParticle* thePrim = 0;
  for (size_t i=0; i<particles.size(); i++) {
    const Data& p = particles[i];	// For convenience below

    thePrim = new G4PrimaryParticle();
    thePrim->SetParticleDefinition(p.pd);
    thePrim->SetMomentumDirection(p.dir);
    thePrim->SetKineticEnergy(p.ekin);
    thePrim->SetWeight(p.wt);
    primaries.push_back(thePrim);

    if (verboseLevel==3) {
      G4cout << i << " : " << p.pd->GetParticleName() << " " << p.ekin/eV
	     << " eV along " << p.dir << " (wt " << p.wt << ")"
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

  // Store position information in summary block
  if (summary) {
    summary->position[0] = pos[0];
    summary->position[1] = pos[1];
    summary->position[2] = pos[2];
    summary->position[3] = time;
  }

  std::vector<G4PrimaryParticle*> primaries;	// Can we make this mutable?
  GetPrimaries(primaries);

  G4double chargeEtot = 0.;		// Cumulative buffers for diagnostics
  G4double phononEtot = 0.;
  G4double bandgap = GetLattice()->GetBandGapEnergy()/2.;

  // Get volume touchable at point and enforce "IsInside()" position
  G4VTouchable* touch = G4CMP::CreateTouchableAtPoint(pos);
  G4ThreeVector newpos = G4CMP::ApplySurfaceClearance(touch, pos);

  // Generate charge carriers in region around track position
  G4bool doCloud = G4CMPConfigManager::CreateChargeCloud();	// Convenience
  if (doCloud) {
    cloud->SetVerboseLevel(verboseLevel);
    cloud->SetTouchable(touch);
    cloud->Generate(nPairsGen, newpos);
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
      G4ThreeVector binpos = chgbin>=0 ? cloud->GetBinCenter(chgbin) : newpos;
      vertex = CreateVertex(event, newpos, time);
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


// Fill a pre-defined set of positions with primaries (no charge cloud)

void G4CMPEnergyPartition::
GetPrimaries(G4Event* event, const std::vector<G4ThreeVector>& pos,
	     G4double time, G4int maxPerVertex) const {
  if (verboseLevel) {
    G4cout << "G4CMPEnergyPartition::GetPrimaries across " << pos.size()
	   << " positions, up to " << maxPerVertex << " per vertex" << G4endl;
  }

  // Use average position in set for summary
  size_t npos = pos.size();
  G4ThreeVector avgpos;
  for (auto const& ip: pos) avgpos += ip;
  avgpos /= npos;

  if (summary) {
    summary->position[0] = avgpos[0];
    summary->position[1] = avgpos[1];
    summary->position[2] = avgpos[2];
    summary->position[3] = time;
  }

  // Create and fill a vertex for each position in set
  size_t tracksPerPos = GetNumberOfTracks() / npos;
  size_t extraTracks = GetNumberOfTracks() - (tracksPerPos * npos);

  std::vector<G4PrimaryParticle*> primaries;	// Can we make this mutable?
  GetPrimaries(primaries);

  size_t iprim = 0;
  for (size_t i=0; i<npos; i++) {
    size_t nprim = tracksPerPos + (i<extraTracks?1:0);
    if (verboseLevel>1)
      G4cout << " adding " << nprim << " primaries to vertex " << i << G4endl;

    G4PrimaryVertex* vtx = 0;

    for (size_t j=0; j<nprim; j++) {
      if (!vtx ||
	  (maxPerVertex>0 && vtx->GetNumberOfParticle()>maxPerVertex)) {
	vtx = CreateVertex(event, pos[i], time);
      }
      vtx->SetPrimary(primaries[iprim++]);
    }	// for (j
  }	// for (i
}


// Create primary vertex at specified location for filling

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

  // If particle type not already set, get it from the current track
  if (summary && summary->PDGcode == 0)
    summary->PDGcode = GetCurrentParticle()->GetPDGEncoding();

  // Store position information in summary block, if not already done
  if (summary && summary->trackID == 0) {
    summary->position[0] = GetCurrentTrack()->GetPosition()[0];
    summary->position[1] = GetCurrentTrack()->GetPosition()[1];
    summary->position[2] = GetCurrentTrack()->GetPosition()[2];
    summary->position[3] = GetCurrentTrack()->GetGlobalTime();
    
    summary->trackID = GetCurrentTrack()->GetTrackID();
    summary->stepID  = GetCurrentTrack()->GetCurrentStepNumber();
  }
  
  // Pre-allocate buffer for secondaries
  secondaries.clear();
  secondaries.reserve(particles.size());

  // Generate charge carriers in region around track position
  G4bool doCloud = G4CMPConfigManager::CreateChargeCloud();	// Convenience
  if (doCloud) {
    cloud->SetVerboseLevel(verboseLevel);
    cloud->SetTouchable(GetCurrentTouchable());
    cloud->Generate(nPairsGen, GetCurrentTrack()->GetPosition());
  }

  if (verboseLevel>1) G4cout << " processing " << particles.size() << G4endl;

  G4Track* theSec = 0;
  G4int ichg = 0;			// Index to deal with charge cloud

  for (size_t i=0; i<particles.size(); i++) {
    const Data& p = particles[i];	// For convenience below

    // Set weights so that generated particles map back to expected true number
    theSec = G4CMP::CreateSecondary(*GetCurrentTrack(), p.pd, p.dir, p.ekin);
    theSec->SetWeight(trkWeight*p.wt);
    secondaries.push_back(theSec);

    // Adjust positions of charges according to generated distribution
    if (doCloud && G4CMP::IsChargeCarrier(theSec))
      theSec->SetPosition(cloud->GetPosition(ichg++));

    if (verboseLevel==3) {
      G4cout << i << " : " << p.pd->GetParticleName() << " " << p.ekin/eV
	     << " eV along " << p.dir << " (wt " << p.wt << ")"
	     << G4endl;
    } else if (verboseLevel>3) {
      G4cout << i << " : ";
      theSec->GetDynamicParticle()->DumpInfo();	// G4Track has no dump function
      G4cout << "   Track Weight = " << theSec->GetWeight() << G4endl;
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


/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPEnergyPartition.hh
/// \brief Definition of the G4CMPEnergyPartition class
///   Functionality to convert energy deposition from Geant4 (both total
///   and non-ionizing) into phonons and charge carrier pairs
///
// $Id$
//
// 20170524  Add constructor and accessor for position argument
// 20170525  Add "rule of five" copy/move operators
// 20170802  Add constructor and accessor for volume argument, particle change
// 20170830  Add function to compute downsampling factors for input energy
// 20170901  Add support for putting primaries directly into event
// 20170925  Add support for distributing charges around position
// 20180424  Need default ctor for Data to support vector::resize()
// 20180425  Add minimum particle generation for downsampling
// 20180827  Add flag to suppress use of downsampling energy scale
// 20190714  Pass particle information through to NuclearRecoil, Lindhard
// 20200218  Support writing DoPartion() internals to event summary data
// 20200222  Add control flag to turn off creating summary data
// 20200805  Add bias across volume to estimate Luke gain downsampling
// 20210328  Split ComputeDownsampling() into individual computation functions
// 20210820  Rename particle count data member for clarity, add counts for
//		after downsampling.  Store weight for each particle in "Data".
// 20220216  Add interface to do partitioning directly from StepAccumulator.
// 20220816  Add generated track counts, for convenience before filling
// 20220816  G4CMP-308 -- Support generating multiple primary positions.

#ifndef G4CMPEnergyPartition_hh
#define G4CMPEnergyPartition_hh 1

#include "globals.hh"
#include "G4CMPProcessUtils.hh"
#include "G4ThreeVector.hh"
#include <vector>

class G4CMPChargeCloud;
class G4CMPPartitionData;
class G4CMPStepAccumulator;
class G4Event;
class G4LatticePhysical;
class G4Material;
class G4ParticleDefinition;
class G4PrimaryParticle;
class G4PrimaryVertex;
class G4Track;
class G4VParticleChange;
class G4VPhysicalVolume;


class G4CMPEnergyPartition : public G4CMPProcessUtils {
public:
  G4CMPEnergyPartition(G4Material* mat=0, G4LatticePhysical* lat=0);
  explicit G4CMPEnergyPartition(const G4VPhysicalVolume* volume);
  explicit G4CMPEnergyPartition(const G4ThreeVector& pos);

  virtual ~G4CMPEnergyPartition();

  // Default copy and move operators
  G4CMPEnergyPartition(const G4CMPEnergyPartition&) = default;
  G4CMPEnergyPartition(G4CMPEnergyPartition&&) = default;
  G4CMPEnergyPartition& operator=(const G4CMPEnergyPartition&) = default;
  G4CMPEnergyPartition& operator=(G4CMPEnergyPartition&&) = default;

  // Set debugging output
  void SetVerboseLevel(G4int vb) { verboseLevel = vb; }

  // Enable or disable summary data collection
  void FillSummary(G4bool fill) { fillSummaryData = fill; }
  G4bool FillingSummary() const { return fillSummaryData; }

  // Toggle whether or not to apply downsampling scale calculations
  void UseDownsampling(G4bool value) { applyDownsampling = value; }
  G4bool UseDownsampling() const { return applyDownsampling; }

  // Placement volume may be used to get material and lattice
  void UseVolume(const G4VPhysicalVolume* volume);

  // Position may be used to get material and lattice from geometry
  void UsePosition(const G4ThreeVector& pos);

  // Material is needed for (Z,A) in Lindhard scaling
  void SetMaterial(G4Material* mat) { material = mat; }

  // Bias voltage may be used to estimate energy from Luke gain
  void SetBiasVoltage(G4double v) { biasVoltage = v; }
  void SetBiasVoltage(const G4ThreeVector& pos);

  // Specify particle type (PDG), total and NIEL energy deposit
  void DoPartition(G4int PDGcode, G4double energy, G4double eNIEL);

  // Specify container with information from one or more G4 steps
  void DoPartition(const G4CMPStepAccumulator* accumulator);

  // Nuclear recoil deposit uses Lindhard scale factor for e/h vs. phonons
  void NuclearRecoil(G4double energy, G4double Z, G4double A);

  // Pure ionization produces no phonons
  void Ionization(G4double energy) { DoPartition(energy, 0.); }

  // Some processes can specify non-ionizing energy directly
  void DoPartition(G4double eIon, G4double eNIEL);

  // Return either primary or secondary particles from partitioning
  void GetPrimaries(std::vector<G4PrimaryParticle*>& primaries) const;

  void GetPrimaries(G4Event* event, const G4ThreeVector& pos, G4double time,
		    G4int maxPerVertex=100000) const;

  void GetPrimaries(G4Event* event, const std::vector<G4ThreeVector>& pos,
		    G4double time, G4int maxPerVertex=100000) const;

  void GetSecondaries(std::vector<G4Track*>& secondaries,
		      G4double trkWeight=1.) const;

  void GetSecondaries(G4VParticleChange* aParticleChange) const;

  // Return number of generated tracks, for convenience before filling vectors
  size_t GetNumberOfPhonons() const { return nPhononsGen; }
  size_t GetNumberOfCharges() const { return nPairsGen*2; }
  size_t GetNumberOfTracks() const { return nPhononsGen + nPairsGen*2; }

  // Assign energy-dependent sampling factors for phonons and charge carriers
  void ComputeDownsampling(G4double eIon, G4double eNIEL);
  void ComputeChargeSampling(G4double eIon);
  void ComputePhononSampling(G4double eNIEL);
  void ComputeLukeSampling(G4double eIon);
  
  // Fraction of total energy deposit in material which goes to e/h pairs
  G4double LindhardScalingFactor(G4double energy, G4double Z=0,
				 G4double A=0) const;

  // Number of e/h pairs to generate including Fano fluctuations
  G4double MeasuredChargePairs(G4double eTrue) const;

protected:
  void GenerateCharges(G4double energy);
  void AddChargePair(G4double ePair, G4double wt);

  void GeneratePhonons(G4double energy);
  void AddPhonon(G4double ePhon, G4double wt);

  G4PrimaryVertex* CreateVertex(G4Event* event, const G4ThreeVector& pos,
				G4double time) const;

  // Create buffer save DoPartition() computations
  G4CMPPartitionData* CreateSummary();

protected:
  G4int verboseLevel;		// Higher numbers give more details
  G4bool fillSummaryData;	// Fill G4CMPPartitionSummary if set

  G4Material* material;		// To get (Z,A) for Lindhard scaling
  G4double biasVoltage;		// Bias across volume for Luke downsampling
  G4double holeFraction;	// Energy from e/h pair taken by hole (50%)
  G4int nParticlesMinimum;	// Minimum production when downsampling
  G4bool applyDownsampling;	// Flag whether to do downsampling calcualtions

  G4CMPChargeCloud* cloud;	// Distribute e/h around central position

  size_t nPairsTrue;		// True number of pairs (no downsampling)
  size_t nPairsGen;		// Number of pairs after downsampling
  G4double chargeEnergyLeft;	// Energy to partition into e/h pairs

  size_t nPhononsTrue;		// True number of phonons (no downsampling)
  size_t nPhononsGen;		// Number of direct phonons after downsampling
  G4double phononEnergyLeft;	// Energy to partition into phonons

  G4CMPPartitionData* summary;	// Summary block, saved to G4HitsCollection

  static const G4ThreeVector origin;
  struct Data {
    G4ParticleDefinition* pd;
    G4ThreeVector dir;
    G4double ekin;
    G4double wt;

    Data() : pd(0), ekin(0.), wt(0.) {;}	// Default ctor for vector::resize()
    Data(G4ParticleDefinition* part, const G4ThreeVector& d, G4double E,
	 G4double w) : pd(part), dir(d), ekin(E), wt(w) {;}
  };
    
  std::vector<Data> particles;	// Combined phonons and charge carriers
};

#endif	/* G4CMPEnergyPartition_hh */

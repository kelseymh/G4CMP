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

#ifndef G4CMPEnergyPartition_hh
#define G4CMPEnergyPartition_hh 1

#include "globals.hh"
#include "G4CMPProcessUtils.hh"
#include "G4ThreeVector.hh"
#include <vector>

class G4CMPChargeCloud;
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

  // Placement volume may be used to get material and lattice
  void UseVolume(const G4VPhysicalVolume* volume);

  // Position may be used to get material and lattice from geometry
  void UsePosition(const G4ThreeVector& pos);

  // Material is needed for (Z,A) in Lindhard scaling
  void SetMaterial(G4Material* mat) { material = mat; }

  // Specify particle type (PDG), total and NIEL energy deposit
  void DoPartition(G4int PDGcode, G4double energy, G4double eNIEL);

  // Nuclear recoil deposit uses Lindhard scale factor for e/h vs. phonons
  void NuclearRecoil(G4double energy) {
    G4double lind = LindhardScalingFactor(energy);
    DoPartition(energy*lind, energy*(1.-lind));
  }

  // Pure ionization produces no phonons
  void Ionization(G4double energy) { DoPartition(energy, 0.); }

  // Some processes can specify non-ionizing energy directly
  void DoPartition(G4double eIon, G4double eNIEL);

  // Return either primary or secondary particles from partitioning
  void GetPrimaries(std::vector<G4PrimaryParticle*>& primaries) const;

  void GetPrimaries(G4Event* event, const G4ThreeVector& pos, G4double time,
		    G4int maxVertex=100000) const;

  void GetSecondaries(std::vector<G4Track*>& secondaries,
		      G4double trkWeight=1.) const;

  void GetSecondaries(G4VParticleChange* aParticleChange) const;

  // Assign energy-dependent sampling factors for phonons and charge carriers
  void ComputeDownsampling(G4double eIon, G4double eNIEL);

  // Fraction of total energy deposit in material which goes to e/h pairs
  G4double LindhardScalingFactor(G4double energy) const;

  // Portion of ionization energy which goes to e/h pairs (Fano factor)
  G4double MeasuredChargeEnergy(G4double eTrue) const;

protected:
  void GenerateCharges(G4double energy);
  void AddChargePair(G4double ePair);

  void GeneratePhonons(G4double energy);
  void AddPhonon(G4double ePhon);

  G4PrimaryVertex* CreateVertex(G4Event* event, const G4ThreeVector& pos,
				G4double time) const;

protected:
  G4int verboseLevel;		// Higher numbers give more details
  G4Material* material;		// To get (Z,A) for Lindhard scaling
  G4double holeFraction;	// Energy from e/h pair taken by hole (50%)
  G4int nParticlesMinimum;	// Minimum production when downsampling

  G4CMPChargeCloud* cloud;	// Distribute e/h around central position
  size_t nCharges;		// Actual (downsampled) number of e+h for cloud

  size_t nPairs;		// True number of pairs (no downsampling)
  G4double chargeEnergyLeft;	// Energy to partition into e/h pairs

  size_t nPhonons;		// True number of phonons (no downsampling)
  G4double phononEnergyLeft;	// Energy to partition into phonons

  static const G4ThreeVector origin;
  struct Data {
    G4ParticleDefinition* pd;
    G4ThreeVector dir;
    G4double ekin;

    Data() : pd(0), ekin(0.) {;}	// Default ctor for vector::resize()
    Data(G4ParticleDefinition* part, const G4ThreeVector& d, G4double E)
      : pd(part), dir(d), ekin(E) {;}
  };
    
  std::vector<Data> particles;	// Combined phonons and charge carriers
};

#endif	/* G4CMPEnergyPartition_hh */

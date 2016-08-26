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

#include "globals.hh"
#include "G4CMPProcessUtils.hh"
#include "G4ThreeVector.hh"
#include <vector>

class G4LatticePhysical;
class G4Material;
class G4ParticleDefinition;
class G4PrimaryParticle;
class G4Track;


class G4CMPEnergyPartition : public G4CMPProcessUtils {
public:
  G4CMPEnergyPartition(G4Material* mat=0, G4LatticePhysical* lat=0);
  virtual ~G4CMPEnergyPartition();

  // Material is needed for (Z,A) in Lindhard scaling
  void SetMaterial(G4Material* mat) { material = mat; }

  // Set debugging output
  void SetVerboseLevel(G4int vb) { verboseLevel = vb; }

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
  void GetSecondaries(std::vector<G4Track*>& secondaries) const;
  
  // Fraction of total energy deposit in material which goes to e/h pairs
  G4double LindhardScalingFactor(G4double energy) const;

  // Portion of ionization energy which goes to e/h pairs (Fano factor)
  G4double MeasuredChargeEnergy(G4double eTrue) const;

protected:
  void GenerateCharges(G4double energy);
  void AddChargePair(G4double ePair);

  void GeneratePhonons(G4double energy);
  void AddPhonon(G4double ePhon);

protected:
  G4Material* material;		// To get (Z,A) for Lindhard scaling
  G4double holeFraction;	// Energy from e/h pair taken by hole (50%)
  G4int verboseLevel;		// Higher numbers give more details

  size_t nPairs;		// Estimated number of pairs to produce
  G4double chargeEnergyLeft;	// Energy to partition into e/h pairs

  size_t nPhonons;		// Estimated number of phonons to produce
  G4double phononEnergyLeft;	// Energy to partition into phonons

  struct Data {
    G4ParticleDefinition* pd;
    G4ThreeVector dir;
    G4double ekin;

    Data(G4ParticleDefinition* part, const G4ThreeVector& d, G4double E)
      : pd(part), dir(d), ekin(E) {;}
  };
    
  std::vector<Data> particles;	// Combined phonons and charge carriers
};

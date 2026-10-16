/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPQPRecombinationRate.cc
/// \brief Compute rate for QP recombining with ambient bath QP into 2Delta phonon
//
// 20260826  G4CMP-662: Splitting off functionality to check rate lookup table
//           for a given input lattice


#ifndef G4CMPQPRecombinationRate_hh
#define G4CMPQPRecombinationRate_hh 1

#include "G4CMPVScatteringRate.hh"
#include <vector>
#include <map>

class G4CMPQPRecombinationRate : public G4CMPVScatteringRate {
					   
public:
  G4CMPQPRecombinationRate() : G4CMPVScatteringRate("qpRecombination") {;}
  virtual ~G4CMPQPRecombinationRate() {;}
  virtual G4double Rate(const G4Track& aTrack) const;
  virtual G4double Rate(const G4Track& aTrack,
                        const G4LatticePhysical * theLat) const;
  virtual void UpdateLookupTable(const G4LatticePhysical * theLat);
  virtual bool CheckLookupTableForLat(const G4LatticePhysical * theLat) const;
  
private:
  
  //Lookup tables for calculated quantities
  std::map<const G4LatticePhysical*,std::vector<std::vector<G4double> > >
  fMap_physicalLattice_NormalizedTauRecombinationVsEnergy;
  std::vector<std::vector<G4double> >
  fCurrentNormalizedTauRecombinationVsEnergy;
  
  //This one doesn't need to be public or protected
  std::vector<std::vector<G4double> >
  ComputeNormalizedTauRecombinationVsEnergy();

  bool CheckToSeeSCParametersSet(const G4LatticePhysical * theLat = 0) const;
  
  //For testing purposes
  void SaveBogoliubovRecombinationTauVsPhononEnergyToLogFile(std::vector<std::vector<G4double> > theFunc);
  
  
};

#endif	/* G4CMPQPRecombinationRate_hh */

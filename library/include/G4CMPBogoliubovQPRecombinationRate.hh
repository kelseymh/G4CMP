/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPBogoliubovQPRecombinationRate.cc
/// \brief Compute rate for QP recombining with ambient bath QP into 2Delta phonon
//

#ifndef G4CMPBogoliubovQPRecombinationRate_hh
#define G4CMPBogoliubovQPRecombinationRate_hh 1

#include "G4CMPVScatteringRate.hh"
#include "G4CMPSCUtils.hh"
#include <vector>
#include <map>

class G4CMPBogoliubovQPRecombinationRate : public G4CMPVScatteringRate {
					   
public:
  G4CMPBogoliubovQPRecombinationRate() : G4CMPVScatteringRate("BogoliubovQPRecombination") {;}
  virtual ~G4CMPBogoliubovQPRecombinationRate() {;}
  
  virtual G4double Rate(const G4Track& aTrack) const;
  virtual void UpdateLookupTable(const G4LatticePhysical * theLat);
  

private:
  
  //Lookup tables for calculated quantities
  std::map<const G4LatticePhysical*,std::vector<std::vector<G4double> > > fMap_physicalLattice_NormalizedTauRecombinationVsEnergy;
  std::vector<std::vector<G4double> > fCurrentNormalizedTauRecombinationVsEnergy;
  
  //This one doesn't need to be public or protected
  std::vector<std::vector<G4double> > ComputeNormalizedTauRecombinationVsEnergy();

  bool CheckToSeeSCParametersSet() const;
  
  //For testing purposes
  void SaveBogoliubovRecombinationTauVsPhononEnergyToLogFile(std::vector<std::vector<G4double> > theFunc);
  
  
};

#endif	/* G4CMPBogoliubovQPRecombinationRate_hh */

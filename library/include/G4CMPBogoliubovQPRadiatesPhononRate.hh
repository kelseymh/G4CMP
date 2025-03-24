/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPBogoliubovQPRadiatesPhononRate.cc
/// \brief Compute rate for QP radiating individual phonons.
//

#ifndef G4CMPBogoliubovQPRadiatesPhononRate_hh
#define G4CMPBogoliubovQPRadiatesPhononRate_hh 1

#include "G4CMPVScatteringRate.hh"
#include "G4CMPSCUtils.hh"
#include <vector>
#include <map>

class G4CMPBogoliubovQPRadiatesPhononRate : public G4CMPVScatteringRate {
					   
public:
  G4CMPBogoliubovQPRadiatesPhononRate() : G4CMPVScatteringRate("BogoliubovQPRadiatesPhonon") {;}
  virtual ~G4CMPBogoliubovQPRadiatesPhononRate() {;}
  
  virtual G4double Rate(const G4Track& aTrack) const;
  virtual void UpdateLookupTable(const G4LatticePhysical * theLat);
  

private:
  
  //Lookup tables for calculated quantities
  std::map<const G4LatticePhysical*,std::vector<std::vector<G4double> > > fMap_physicalLattice_NormalizedTauQPRadiatesPhononVsEnergy;
  std::vector<std::vector<G4double> > fCurrentNormalizedTauQPRadiatesPhononVsEnergy;
  
  //This one doesn't need to be public or protected
  std::vector<std::vector<G4double> > ComputeNormalizedTauQPRadiatesPhononVsEnergy();
  
  bool CheckToSeeSCParametersSet() const;
  
  //For testing purposes
  void SaveBogoliubovQPRadiatesPhononTauVsPhononEnergyToLogFile(std::vector<std::vector<G4double> > theFunc);
  
  
};

#endif	/* G4CMPBogoliubovQPRadiatesPhononRate_hh */

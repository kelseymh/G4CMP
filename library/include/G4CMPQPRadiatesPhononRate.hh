/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPQPRadiatesPhononRate.cc
/// \brief Compute rate for QP radiating individual phonons.
//

#ifndef G4CMPQPRadiatesPhononRate_hh
#define G4CMPQPRadiatesPhononRate_hh 1

#include "G4CMPVScatteringRate.hh"
#include "G4CMPSCUtils.hh"
#include <vector>
#include <map>

class G4CMPQPRadiatesPhononRate : public G4CMPVScatteringRate {
					   
public:
  G4CMPQPRadiatesPhononRate() : G4CMPVScatteringRate("QPRadiatesPhonon") {;}
  virtual ~G4CMPQPRadiatesPhononRate() {;}
  
  virtual G4double Rate(const G4Track& aTrack) const;
  virtual void UpdateLookupTable(const G4LatticePhysical * theLat);
  

private:
  
  //Lookup tables for calculated quantities
  std::map<const G4LatticePhysical*,std::vector<std::vector<G4double> > >
  fMap_physicalLattice_NormalizedTauQPRadiatesPhononVsEnergy;
  std::vector<std::vector<G4double> >
  fCurrentNormalizedTauQPRadiatesPhononVsEnergy;
  
  //This one doesn't need to be public or protected
  std::vector<std::vector<G4double> >
  ComputeNormalizedTauQPRadiatesPhononVsEnergy();
  
  bool CheckToSeeSCParametersSet() const;
  
  //For testing purposes
  void SaveQPRadiatesPhononTauVsPhononEnergyToLogFile(std::vector<std::vector<G4double> > theFunc);
  
  
};

#endif	/* G4CMPQPRadiatesPhononRate_hh */

/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPPairBreakingRate.hh
/// \brief Compute rate for phonon breaking CPs into 2 QPs.
///
/// This class computes and dictates the rate for phonons breaking
/// Cooper pairs into quasiparticles in superconductors. Paired with
/// G4CMPSCPairbreakingProcess.
//
// 20250922  G4CMP-219 -- First addition to this history (done at time
//                        of merge to develop)

#ifndef G4CMPSCPairBreakingRate_hh
#define G4CMPSCPairBreakingRate_hh 1

#include "G4CMPVScatteringRate.hh"
#include "G4CMPSCUtils.hh"

class G4CMPSCPairBreakingRate : public G4CMPVScatteringRate {  
public:
  G4CMPSCPairBreakingRate() : G4CMPVScatteringRate("scPairBreaking") {;}
  virtual ~G4CMPSCPairBreakingRate() {;}

  virtual G4double Rate(const G4Track& aTrack) const;
  virtual void UpdateLookupTable(const G4LatticePhysical* theLat);

private:

  //Lookup tables for calculated quantities
  std::map<const G4LatticePhysical*,std::vector<std::vector<G4double> > >
  fMap_physicalLattice_NormalizedTauPairBreakingVsEnergy;
  
  std::vector<std::vector<G4double> >
  fCurrentNormalizedTauPairBreakingVsEnergy;
  
  //This one doesn't need to be public or protected
  std::vector<std::vector<G4double> >
  ComputeNormalizedTauPairBreakingVsEnergy();

  bool CheckToSeeSCParametersSet() const;  
  void SavePairBreakingRateVsPhononEnergyToLogFile(std::vector<std::vector<G4double> > theFunc);
};

#endif	/* G4CMPSCPairBreakingRate_hh */

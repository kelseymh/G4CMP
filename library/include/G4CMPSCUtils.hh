/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPSCUtils.hh
/// \brief Definition of the G4CMPSCUtils class
///   Provides useful general functions for use by KaplanSCPhonon and BogoliubovQP processes
///

#ifndef G4CMPSCUtils_h
#define G4CMPSCUtils_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include <map>
#include <utility>
#include <vector>

class G4LatticePhysical;


class G4CMPSCUtils {
public:
  G4CMPSCUtils();
  virtual ~G4CMPSCUtils();

  G4CMPSCUtils(const G4CMPSCUtils&) = default;
  G4CMPSCUtils(G4CMPSCUtils&&) = default;
  G4CMPSCUtils& operator=(const G4CMPSCUtils&) = default;
  G4CMPSCUtils& operator=(G4CMPSCUtils&&) = default;
  
  virtual void SetVerboseLevel(G4int vb) { scuVerboseLevel = vb; }






protected:

  //Processing: called by process and rate classes, so should be at least protected
  void SetCurrentSCInfoToNull();  
  void LoadLatticeInfoIntoSCUtils(const G4LatticePhysical * theLat);
  
  //Since now rate classes and processes are actually the ones doing the calculation, we need to make these protected instead of private
  G4double FermiFactor(G4double energy, G4double temperature);
  G4double BoseFactor(G4double energy, G4double temperature);
  G4double GetTauAsAFunctionOfEnergy( const std::vector<std::vector<G4double> > & tauVsPhononEnergy, G4String particleInQuestion, G4double energy ) const;
  G4double ComputeTestGapEnergyAtNonzeroT(double Teff, double Tcrit, double gap0Energy) const;
  
private:

  //These run fully internally -- can be private
  G4double ComputeCurrentGapEnergyAtNonzeroT();


  
  G4int scuVerboseLevel;			// For local use; name avoids collisions

protected:

  //Want all of the SC info to be accessible to the rate classes, which inherit from the G4CMPSCUtils,
  //so this should be protected. Actually, should migrate away from this eventually and do get functions
  //so that rate classes can't modify these... REL. That's the right way to do this.
  
  //-------------------------------------------
  //SC parameters needed for computations
  //Stored/input parameters
  G4double fGap0Energy;                       //Gap energy at T=0
  G4double fTau0_qp;
  G4double fTau0_ph;
  G4double fSoundSpeed;
  G4double fTcrit;
  G4double fTeff;
  G4double fDn;
    
  //Computed parameters
  G4double fGapEnergy;

  //-------------------------------------------
  //Parameters for defining the various lookup tables needed by this class. 
  enum{ fPhononEnergyBins=2000, fQPEnergyBins=2000, fGapEnergyTempDependenceBins=53 }; //Originally 1000,1000,53
  G4double fMinPhononEnergyDivGap;
  G4double fMaxPhononEnergyDivGap;
  G4double fMinQPEnergyDivGap;
  G4double fMaxQPEnergyDivGap;
  
  //-------------------------------------------
  //Maps and lookup tables
  G4double fGapEnergyTempDependence[fGapEnergyTempDependenceBins][2];
  

};

#endif	/* G4CMPSCUtils_hh */

/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPEmpiricalLindhardNIEL.hh
/// \brief NIEL calculation from Ionization yield measurement 
/// in a germanium CDMSlite detector using photo-neutron sources.
///
/// Computation of NIEL according to Lindhard (Lewin Smith) yield with 
/// a parameter k that depends on the energy.
///
/// Paper DOI: https://doi.org/10.1103/PhysRevD.105.122002
//
// 20250212  David Sadek

#ifndef G4CMPEmpiricalLindhardNIEL_hh
#define G4CMPEmpiricalLindhardNIEL_hh 1

#include "G4VNIELPartition.hh"
#include "G4CMPLewinSmithNIEL.hh"
#include "G4SystemOfUnits.hh"


class G4CMPEmpiricalLindhardNIEL : public G4CMPLewinSmithNIEL {
public:
    // Constructor with default values
  G4CMPEmpiricalLindhardNIEL(G4double klow_ = 0.040, 
                             G4double khigh_ = 0.142, 
                             G4double Elow_ = 0.39 * keV, 
                             G4double Ehigh_ = 7.0 * keV, 
                             bool useEnergyDependentK_ = true, 
                             G4double kFixed_ = 0.158)
    : klow(klow_), khigh(khigh_), Elow(Elow_), Ehigh(Ehigh_), 
      useEnergyDependentK(useEnergyDependentK_), kFixed(kFixed_) {;}
  virtual ~G4CMPEmpiricalLindhardNIEL() {;}
  
  // return the fraction of the specified energy which will be deposited as NIEL
  // if an incoming particle is stopped in the specified material,
  // energy in native G4 energy units.
  //
  // The parameter k is a function of energy

  // Setters to modify parameters if needed
  void SetKlow(G4double kl) { klow = kl; }
  void SetKhigh(G4double kh) { khigh = kh; }
  void SetElow(G4double el) { Elow = el; }
  void SetEhigh(G4double eh) { Ehigh = eh; }
  void SetFixedK(G4double fixedK) { kFixed = fixedK; }
  void UseEnergyDependentK(bool use) { useEnergyDependentK = use; }

  virtual G4double 
  PartitionNIEL(G4double energy, const G4Material *material, G4double Zin=0.,
		G4double Ain=0.) const override;

private:

  G4double klow;      // Lower value of k
  G4double khigh;     // High value of k
  G4double Elow;      // Lowest energy value 
  G4double Ehigh;     // Highest energy value
  bool useEnergyDependentK;         // Key to determine whether to use energy dependent k
  G4double kFixed;    // Value of k if not energy dependent nor passed
  mutable bool firstCall = true; 	// Print warning messages only once
  mutable bool firstCall_E = true; 	// Print warning messages only once


};

#endif	/* G4CMPEmpiricalLindhardNIEL_hh */

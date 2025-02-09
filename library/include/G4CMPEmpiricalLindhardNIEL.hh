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
// 20250208  David Sadek

#ifndef G4CMPEmpiricalLindhardNIEL_hh
#define G4CMPEmpiricalLindhardNIEL_hh 1

#include "G4VNIELPartition.hh"
#include "LewinSmithNIEL.hh"


class G4CMPEmpiricalLindhardNIEL : public G4CMPLewinSmithNIEL {
public:
  G4CMPEmpiricalLindhardNIEL() {;}
  virtual ~G4CMPEmpiricalLindhardNIEL() {;}
  
  // return the fraction of the specified energy which will be deposited as NIEL
  // if an incoming particle is stopped in the specified material,
  // energy in native G4 energy units.
  //
  // The parameter k is a function of energy

  virtual G4double 
  PartitionNIEL(G4double energy, const G4Material *material = "G4_Ge", G4double Zin=0.,
		G4double Ain=0., G4double klow = 0.040, G4double khigh = 0.142,                            
  		G4double Elow = 0.39 * keV, G4double Ehigh = 7.0 * keV) const;

private:

  G4double Klow;      // Lower value of k
  G4double Khigh;     // High value of k
  G4double Elow;      // Lowest energy value 
  G4double Ehigh;     // Highest energy value
  mutable bool firstCall = true; 	// Print warning messages only once
  mutable bool firstCall_E = true; 	// Print warning messages only once


};

#endif	/* G4CMPEmpiricalLindhardNIEL_hh */

/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPSarkisNIEL.hh
/// \brief Non-ionizing energy loss calculation from Sarkis 2023.
///
/// Link to the paper:https://journals.aps.org/pra/abstract/10.1103/PhysRevA.107.062811.
///
/// This ionization model was obtained from the Sarkis paper referenced above 
/// for Silicon ONLY so it doos not have (Z,A) dependence. The code will check
/// the effective Z of the input material. If the effZ is within
/// +/-1 of Silicon Z and A, the Sarkis model will be used, else,
/// Lindhard(LewinSmith) model will be used for NIEL calculations.
///
/// The model is obtained in the range of ~50 eV to 3 MeV above which
/// Lindahard (LewinSmith) model will be used for NIEL calculations

// 20230721  David Sadek - University of Florida (david.sadek@ufl.edu)
// 20250102  M. Kelsey -- Move constructor implementation to .cc file.
// 20250211  D. Sadek Fix the energy scale and interpolation.

#ifndef G4CMPSarkisNIEL_hh
#define G4CMPSarkisNIEL_hh 1

#include "G4CMPLewinSmithNIEL.hh"
#include "G4PhysicsFreeVector.hh"


class G4CMPSarkisNIEL : public G4CMPLewinSmithNIEL {
public:
  G4CMPSarkisNIEL();
  virtual ~G4CMPSarkisNIEL() {;}
  
  // return the fraction of the specified energy which will be deposited as NIEL
  // if an incoming particle with z1, a1 is stopped in the specified material
  // a1 is in atomic mass units, energy in native G4 energy units.
  //
  virtual G4double 
  PartitionNIEL(G4double energy, const G4Material * material, G4double Zin = 0.,
		G4double Ain = 0.) const override;
    
private:
  const G4double SiZ = 14.0;
  const G4double SiA = 28.09;
  
  mutable G4PhysicsFreeVector lVector;
  mutable std::size_t idx = 0;		// Dummy argument to G4PLV::Value()
  mutable bool firstCall = true; 	// Print warning messages only once
    
};

#endif	/* G4CMPSarkisNIEL_hh */

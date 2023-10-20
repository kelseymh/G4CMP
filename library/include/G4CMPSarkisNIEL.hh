/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPSarkisNIEL.hh
/// \brief Non-ionizing energy loss calculation from Sarkis 2023.
///
/// Link to the paper:https://journals.aps.org/pra/abstract/10.1103/PhysRevA.107.062811.
//

// 20230721  David Sadek - University of Florida (david.sadek@ufl.edu)

// This ionization model was obtained from the Sarkis paper referenced above for Silicon ONLY so it deos not have (Z,A) dependence. The code will check the effective Z and A of the input material, if the effZ and effA are within +/-1 of Silicon Z and A, the Sarkis model will be used, else, Lindhard(LewinSmith) model will be used for NIEL calculations.


// The model is obtained in the range of ~50 eV to 3 MeV above which Lindahard(LewinSmith) model will be used for NIEL calculations

#ifndef G4CMPSarkisNIEL_hh
#define G4CMPSarkisNIEL_hh 1

#include "G4CMPLewinSmithNIEL.hh"
#include "G4CMPConfigManager.hh"
#include "G4PhysicsLinearVector.hh"

class G4CMPSarkisNIEL : public G4CMPLewinSmithNIEL {
public:
  G4CMPSarkisNIEL() : fPath(G4CMPConfigManager::GetLatticeDir() + "/NIEL/SarkisSiIonYield.txt") {;} //inputFile.open(fPath)
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
    G4String fPath;                     // full path to the data file
    //std::ifstream inputFile;
    mutable G4PhysicsLinearVector lVector;
    mutable std::size_t idx = 0;
    mutable bool firstCall = true; // A static variable to be used to print a warning message only once if the material passed is not Silicon
    
};

#endif	/* G4CMPSarkisNIEL_hh */

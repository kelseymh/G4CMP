/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPImpactTUNLNIEL.hh
/// \brief Non-ionizing energy loss calculation from IMPACT@TUNL 2023.
///
/// Computation of NIEL using the empirical model extracted from the IMPACT@TUNL ionization yield measurements.  Link to the paper:https://arxiv.org/abs/2303.02196.
//

// 20230721  David Sadek  University of Florida (david.sadek@ufl.edu)

// This ionization model was obtained from the ionization yield measurements in Silicon ONLY and it deos not have (Z,A) dependence. The code will check the effective Z and A of the input material, the effZ and effA are within +/-1 of Silicon Z and A, the Impact model will be used, else, Lindhard(LewinSmith) model will be used for NIEL calculations. 

// The model is obtained in the range of 100 eV to 10 keV. Above 10 keV, Lindhard(LewinSmith) model will be used.


#include "globals.hh"
#include "G4CMPImpactTunlNIEL.hh"
#include "G4ExceptionSeverity.hh"
#include "G4Material.hh"
#include "G4Pow.hh"
#include "G4SystemOfUnits.hh"
#include "G4LewinSmithNIEL.hh"
#include <algorithm>
#include <cmath>


// G4CMPImpactTunlNIEL::G4CMPImpactTunlNIEL(){
//     G4Exception("G4CMPImpactTunlNIEL", "G4CMP1003", JustWarning,
//                 "This model is obtained in the range of 100 eV to 10 keV. Above 10 keV, Lindhard model will be used.");
// }

G4CMPImpactTunlNIEL::G4CMPImpactTunlNIEL() : B(0.261), Y10(0.302), SiZ(14.0), SiA(28.09) {}

G4double G4CMPImpactTunlNIEL::
PartitionNIEL(G4double energy, const G4Material *material,G4double Zin=0.,
		G4double Ain=0.) const {
  if (!material) {
    G4Exception("G4CMPImpactTunlNIEL", "G4CMP1000", FatalErrorInArgument,
		  "No material passed to partition function");
    return 1.;		// Won't get here after FATAL error
  }
  // Get effective Z,A of the material
  const G4double Z = GetEffectiveZ(material);
  const G4double A = GetEffectiveA(material) / (g / mole);
    
  // Check if the material is silicon or similar (within +-1 of Z and A of silicon)
    
  if (G4doubleabs(Z - SiZ) <= 1.0 && G4doubleabs(A - SiA) <= 1.0) {
      G4Exception("G4CMPImpactTunlNIEL", "G4CMP1005", JustWarning,
                  "IMPACT@TUNL model is obtained in the range of 100 eV to 10 keV. Above 10 keV, Lindhard model will be used.");
    
    // IMPACT@TUNL model below 10 keV
      
    if (energy <= 10 * keV) {
        G4double EB = g4pow->powA(energy/(10*keV),B);
        return (Y10*EB); 
    } else {
        return G4CMPLewinSmithNIEL:PartitionNIEL(energy, material, Zin, Ain);
    }
  } else {
      G4Exception("G4CMPImpactTunlNIEL", "G4CMP1004", JustWarning, "The input material is not Silicon. The Lindhard model will be used for NIEL calculation.");
      return G4CMPLewinSmithNIEL::PartitionNIEL(energy, material, Zin, Ain);
  }
}

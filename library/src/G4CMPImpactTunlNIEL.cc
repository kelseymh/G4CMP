/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPImpactTUNLNIEL.hh
/// \brief Non-ionizing energy loss calculation from IMPACT@TUNL 2023.
///
/// Computation of NIEL using the empirical model extracted from the ionization yield measurements according to IMPACT@TUNL 2023.  Link to the paper:https://arxiv.org/abs/2303.02196.
//
// $Id$
//
// 20190711  Michael Kelsey

// potential problem1. This ionization model was obtained from the ionization yield measurements in Silicon ONLY.

// potential problem2. The model is obtained in the range of 100 eV to 10 keV. What about energies>10 keV? My first choice would be Lindhard since it works very well above 10 keV but we could give the option to the user to decide which model to use above 10 keV


#include "globals.hh"
#include "G4CMPImpactTunlNIEL.hh"
#include "G4ExceptionSeverity.hh"
#include "G4Material.hh"
#include "G4Pow.hh"
#include "G4SystemOfUnits.hh"
#include "G4LewinSmithNIEL.hh"
#include <algorithm>
#include <cmath>


// return the fraction of the specified energy which will be deposited as NIEL
// if an incoming particle with z1, a1 is stopped in the specified material
// a1 is in atomic mass units, energy in native G4 energy units.
//

G4double G4CMPImpactTunlNIEL::
PartitionNIEL(G4double energy, const G4Material *material,G4double Zin=0.,
		G4double Ain=0.) const {
  if (!material) {
    G4Exception("G4CMPImpactTunlNIEL", "G4CMP1000", FatalErrorInArgument,
		  "No material passed to partition function");
    return 1.;		// Won't get here after FATAL error
    // This could deal with the potential error 1 to exclude any material that's not silicon. Need to discuss how.
  } 
  
  G4Pow* g4pow = G4Pow::GetInstance();		// Fast, tabulated x^(2/3) etc.
    
  // IMPACT@TUNL model below 10 keV 
    
  if (energy<10*keV){
      const G4double B=0.261;  // this is the central value for B=0.261
      const G4double Y10=0.302;
      G4double EB=g4pow->powZ(energy/(10*keV),B);
      return (Y10*EB);
  } else {
      G4CMPLewinSmithNIEL LS;
      return LS.PartitionNIEL(energy, material);
//     // Effective (weighted average) properties of material
    
//       const G4double Z=GetEffectiveZ(material), A=GetEffectiveA(material)/(g/mole);

//     // From Lewin and Smith, 1996
//     // NOTE:  Avoiding use of std::pow for Z^(-7/3) and Z^(2/3) below
//       G4double z23 = g4pow->Z23(Z);
//       G4double epsilon = 11.5/keV * energy / (Z*z23*z23);	// * Z^(-7/3)
//       G4double k = 0.133 * z23 / std::sqrt(A);
//       G4double h = (0.7*g4pow->powA(epsilon,0.6) + 3.*g4pow->powA(epsilon,0.15)
//             + epsilon);

//       return (k*h / (1.+k*h));			// Returns EM fraction
  }
}

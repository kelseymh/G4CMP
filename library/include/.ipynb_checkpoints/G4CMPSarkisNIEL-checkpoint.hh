/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPSarkisNIEL.hh
/// \brief Non-ionizing energy loss calculation from Sarkis 2023.
///
/// Link to the paper:https://journals.aps.org/pra/abstract/10.1103/PhysRevA.107.062811.
//
// $Id: 7dd025e9b6caeb87b0113dec6200bc2c05f9a237 $
//
// 20230721  David Sadek - University of Florida (david.sadek@ufl.edu)


#include "globals.hh"
#include "G4CMPSarkisNIEL.hh"
#include "G4ExceptionSeverity.hh"
#include "G4Material.hh"
#include "G4Pow.hh"
#include "G4SystemOfUnits.hh"
#include "G4InterpolationManager.hh"
#include "G4DataInterpolation.hh"
#include <algorithm>
#include <cmath>


// return the fraction of the specified energy which will be deposited as NIEL
// if an incoming particle with z1, a1 is stopped in the specified material
// a1 is in atomic mass units, energy in native G4 energy units.
//

G4double G4CMPSarkisNIEL::
PartitionNIEL(G4double energy, const G4Material *material,G4double Zin=0.,
		G4double Ain=0.) const {
  if (!material) {
    G4Exception("G4CMPSarkisNIEL", "G4CMP1000", FatalErrorInArgument,
		  "No material passed to partition function");
    return 1.;		// Won't get here after FATAL error
    // This could deal with the potential error 1 to exclude any material that's not silicon. Need to discuss how.
  } 
  
  G4Pow* g4pow = G4Pow::GetInstance();		// Fast, tabulated x^(2/3) etc.
    
  // Sarkis model 
    
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

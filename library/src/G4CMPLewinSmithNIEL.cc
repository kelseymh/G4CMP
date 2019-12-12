/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPLewinSmithNIEL.cc
/// \brief Non-ionizing energy loss calculation from Lewin & Smith 1996.
///
/// Computation of NIEL according to Lewin & Smith 1996, depending only
/// upon material properties.  Used by G4CMPEnergyPartition.
//
// $Id$
//
// 20190711  Michael Kelsey
// 20191211  BUG FIX: Use base class function to get Z,A of compound materials

#include "globals.hh"
#include "G4CMPLewinSmithNIEL.hh"
#include "G4ExceptionSeverity.hh"
#include "G4Material.hh"
#include "G4Pow.hh"
#include "G4SystemOfUnits.hh"
#include <algorithm>
#include <cmath>


// return the fraction of the specified energy which will be deposited as NIEL
// if an incoming particle with z1, a1 is stopped in the specified material
// a1 is in atomic mass units, energy in native G4 energy units.
//
// NOTE: Projectile properties are ignored in Lewin & Smith calculation

G4double G4CMPLewinSmithNIEL::
PartitionNIEL(G4double energy, const G4Material *material, G4double /*Zin*/,
	      G4double /*Ain*/) const {
  if (!material) {
    G4Exception("G4CMPLewinSmithNIEL", "G4CMP1000", FatalErrorInArgument,
		  "No material passed to partition function");
    return 1.;		// Won't get here after FATAL error
  }

  G4Pow* g4pow = G4Pow::GetInstance();		// Fast, tabulated x^(2/3) etc.

  // Effective (weighted average) properties of material
  const G4double Z=GetEffectiveZ(material), A=GetEffectiveA(material)/(g/mole);

  // From Lewin and Smith, 1996
  // NOTE:  Avoiding use of std::pow for Z^(-7/3) and Z^(2/3) below
  G4double z23 = g4pow->Z23(Z);
  G4double epsilon = 0.0115 * energy / (Z*z23*z23);	// * Z^(-7/3)
  G4double k = 0.133 * z23 / std::sqrt(A);
  G4double h = (0.7*std::pow(epsilon,0.6) + 3.*std::pow(epsilon,0.15)
		+ epsilon);

  return (k*h / (1.+k*h));
}

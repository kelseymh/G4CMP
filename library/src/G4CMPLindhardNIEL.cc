/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPLindhardNIEL.cc
/// \brief Non-ionizing energy loss calculation from Lewin & Smith 1996.
///
/// Computation of NIEL according to Lindhard & Robinson, depending upon
/// both the material properties and the projectile mass and charge.
/// May be registered by client code into G4CMPEnergyPartition.
//
// $Id$
//
// 20190711  Michael Kelsey

#include "globals.hh"
#include "G4CMPLindhardNIEL.hh"
#include "G4ExceptionSeverity.hh"
#include "G4Material.hh"
#include "G4Pow.hh"
#include "G4SystemOfUnits.hh"
#include <algorithm>
#include <cmath>


// return the fraction of the specified energy which will be deposited as NIEL
// if an incoming particle with z1, a1 is stopped in the specified material
// a1 is in atomic mass units, energy in native G4 energy units.

G4double G4CMPLindhardNIEL::
PartitionNIEL(G4double energy, const G4Material *material, G4double z1,
	      G4double a1) const {
  if (!material) {
    G4Exception("G4CMPLindhardNIEL", "G4CMP1000", FatalErrorInArgument,
		  "No material passed to partition function");
    return 1.;		// Won't get here after FATAL error
  }

  G4Pow* g4pow = G4Pow::GetInstance();		// Fast, tabulated x^(2/3) etc.

  // Compute partition based on most abundant element of material
  size_t nMatElements = material->GetNumberOfElements();
        
  const G4double* atomDensities = material->GetVecNbOfAtomsPerVolume();
  size_t maxindex = (std::max_element(atomDensities, atomDensities+nMatElements)
		     - atomDensities);

  const G4Element *element = material->GetElement(maxindex);
  G4double z2 = element->GetZ();
  G4double a2 = element->GetA()/(g/mole);
        
  G4double zpow = g4pow->Z23(z1) + g4pow->Z23(z2);
  G4double asum = a1+a2;
        
  G4double el = 30.724*z1*z2*std::sqrt(zpow)*asum/a2;
  G4double fl = (0.0793*g4pow->Z23(z1)*std::sqrt(z2*asum*asum*asum/(a1*a1*a1*a2))
		 / std::pow(zpow, 0.75));
  G4double eps = (energy/eV) / el;

  G4double denom = 1. + fl*(eps + 3.4008*std::pow(eps, 0.16667)
			    + 0.40244*std::pow(eps, 0.75));

  return 1.0/denom;
}

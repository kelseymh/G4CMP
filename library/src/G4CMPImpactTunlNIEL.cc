/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPImpactTUNLNIEL.hh
/// \brief Non-ionizing energy loss calculation from IMPACT@TUNL 2023.
///
/// Computation of NIEL using the empirical model extracted from the
/// IMPACT@TUNL ionization yield measurements.  
/// Link to the paper:https://arxiv.org/abs/2303.02196.
///
/// This ionization model was obtained from the ionization yield measurements
/// in Silicon ONLY and it does not have (Z,A) dependence. The code will
/// check the effective Z of the input material. If the effZ is within
/// +/-1 of Silicon Z, the Impact model will be used, else,
/// Lindhard (LewinSmith) model will be used for NIEL calculations.
///
/// The model is obtained in the range of 100 eV to 10 keV.
/// Above 10 keV, Lindhard (LewinSmith) model will be used.
//
// 20230721  David Sadek - University of Florida (david.sadek@ufl.edu)
// 20240416  S. Zatschler -- Remove unused const A
// 20240417  M. Kelsey -- Add check using SiA and A
// 20250102  M. Kelsey -- Fix placement of warning messages.

#include "globals.hh"
#include "G4CMPImpactTunlNIEL.hh"
#include "G4ExceptionSeverity.hh"
#include "G4Material.hh"
#include "G4Pow.hh"
#include "G4SystemOfUnits.hh"
#include <algorithm>
#include <cmath>


G4double G4CMPImpactTunlNIEL::
PartitionNIEL(G4double energy, const G4Material *material,G4double Zin,
		G4double Ain) const {
  if (!material) {
    G4Exception("G4CMPImpactTunlNIEL", "G4CMP1000", FatalErrorInArgument,
		  "No material passed to partition function");
    return 1.;		// Won't get here after FATAL error
  }

  // Get effective Z of the material
  const G4double Z = GetEffectiveZ(material);
  const G4double A = GetEffectiveA(material);

  // Check if the material is silicon or similar
  if (std::abs(Z-SiZ) < 0.5 || std::abs(A/(g/mole)-SiA) < 0.5) {
    // IMPACT@TUNL model below 10 keV, a simple power law in energy
    G4Pow* g4pow = G4Pow::GetInstance(); 
    if (energy <= 10.*keV) {
      G4double EB = g4pow->powA(energy/(10*keV),B);
      return (Y10*EB); 
    } else {
      if (firstCall) {
	G4Exception("G4CMPImpactTunlNIEL", "G4CMP1005", JustWarning,
		    "IMPACT@TUNL model is obtained in the range of 100 eV to 10 keV.\nAbove 10 keV, LewinSmith model will be used.");
	firstCall = false;
      }
      return G4CMPLewinSmithNIEL::PartitionNIEL(energy, material, Zin, Ain);
    }
  } else {
    if (firstCall) {
      G4cerr << "Z " << Z << " vs SiZ " << SiZ << ", A " << A/(g/mole)
	     << " vs " << " SiA " << SiA << G4endl;

      G4Exception("G4CMPImpactTunlNIEL", "G4CMP1004", JustWarning,
		  "The input material is not Silicon.\nThe LewinSmith model will be used for NIEL calculation.");
      firstCall = false;
    }
    return G4CMPLewinSmithNIEL::PartitionNIEL(energy, material, Zin, Ain);
  }
}

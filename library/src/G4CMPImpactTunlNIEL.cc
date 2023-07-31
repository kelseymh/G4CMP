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

// This ionization model was obtained from the ionization yield measurements in Silicon ONLY and it deos not have (Z,A) dependence. It may not give the correct values for other materials.

// The model is obtained in the range of 100 eV to 10 keV. Above 10 keV, Lindhard model will be used.

// return the fraction of the specified energy which will be deposited as NIEL
// if an incoming particle with z1, a1 is stopped in the specified material
// a1 is in atomic mass units, energy in native G4 energy units.
//

#include "globals.hh"
#include "G4CMPImpactTunlNIEL.hh"
#include "G4ExceptionSeverity.hh"
#include "G4Material.hh"
#include "G4Pow.hh"
#include "G4SystemOfUnits.hh"
#include "G4LewinSmithNIEL.hh"
#include <algorithm>
#include <cmath>


G4CMPImpactTunlNIEL::G4CMPImpactTunlNIEL(){
    G4Exception("G4CMPImpactTunlNIEL", "G4CMP1003", JustWarning,
                "The IMPACT@TUNL ionization yield model was obtained embirically ONLY for Silicon and may not give the correct values for other materials.");
}


G4double G4CMPImpactTunlNIEL::
PartitionNIEL(G4double energy, const G4Material *material,G4double Zin=0.,
		G4double Ain=0.) const {
  if (!material) {
    G4Exception("G4CMPImpactTunlNIEL", "G4CMP1000", FatalErrorInArgument,
		  "No material passed to partition function");
    return 1.;		// Won't get here after FATAL error

  } 
  
  G4Pow* g4pow = G4Pow::GetInstance();		// Fast, tabulated x^(2/3) etc.
    
  // IMPACT@TUNL model below 10 keV 
    
  if (energy<10*keV){
      const G4double B=0.261;  // this is the central value for B=0.261
      const G4double Y10=0.302;
      G4double EB=g4pow->powZ(energy/(10*keV),B);
      return (Y10*EB);
  } 
    // Lindhardh model above 10 keV
  else {
      G4CMPLewinSmithNIEL LS;
      return LS.PartitionNIEL(energy, material);
  }
}

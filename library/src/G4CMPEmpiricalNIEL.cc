/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPEmpiricalNIEL.hh
/// \brief NIEL calculation with energy-dependent k. The function k(E) is
/// taken from the Ionization yield measurement in Ge CDMSlite detector
/// using photo-neutron sources.

/// If energy dependent k is not used, NIEL calculation will be according 
/// to Lindhard (Lewin Smith) yield with k value set by the user.
///
/// Paper DOI: https://doi.org/10.1103/PhysRevD.105.122002
//
// 20250212  David Sade

#include "globals.hh"
#include "G4CMPEmpiricalNIEL.hh"
#include "G4ExceptionSeverity.hh"
#include "G4CMPConfigManager.hh"
#include "G4Material.hh"
#include "G4Pow.hh"
#include "G4SystemOfUnits.hh"
#include <algorithm>
#include <cmath>
#include <cfloat>


// return the fraction of the specified energy which will be deposited as NIEL
// if an incoming particle is stopped in the specified material,
// energy in native G4 energy units.
//
// The parameter k is a function of energy


G4CMPEmpiricalNIEL::G4CMPEmpiricalNIEL() 
    : Empklow(DBL_MIN), Empkhigh(DBL_MIN), 
      EmpElow(DBL_MIN), EmpEhigh(DBL_MIN), 
      EmpkFixed(DBL_MIN), EmpEDepK(true) {

    G4CMPConfigManager* config = G4CMPConfigManager::Instance();

    // Fetch values from G4CMPConfigManager or set defaults
    Empklow = (config->GetEmpklow() != DBL_MIN) ? config->GetEmpklow() : 0.040;
    Empkhigh = (config->GetEmpkhigh() != DBL_MIN) ? config->GetEmpkhigh() : 0.142;
    EmpElow = (config->GetEmpElow() != DBL_MIN) ? config->GetEmpElow() : 0.39 * keV;
    EmpEhigh = (config->GetEmpEhigh() != DBL_MAX) ? config->GetEmpEhigh() : 7.0 * keV;
    EmpkFixed = (config->GetEmpkFixed() != DBL_MIN) ? config->GetEmpkFixed() : 0.158;
    EmpEDepK = config->GetEmpEDepK();
}

G4double G4CMPEmpiricalNIEL::
PartitionNIEL(G4double energy, const G4Material *material, G4double Zin,
		G4double Ain) const {

  if (!material) {
    G4Exception("G4CMPEmpiricalNIEL", "G4CMP1000", FatalErrorInArgument,
		  "No material passed to partition function");
    return 1.;		// Won't get here after FATAL error
  }

  if (useLewinSmith) {
      return G4CMPLewinSmithNIEL::PartitionNIEL(energy, material, Zin, Ain);
    
  if (EmpEDepK && ((energy < EmpElow) || (energy > EmpEhigh))) {
      if (firstCall_E) {
          G4Exception("G4CMPEmpiricalNIEL", "G4CMP1002", JustWarning,
		    "Energy is out of bounds. Lewin-Smith model is used.");
         firstCall_E = false;
      }
      return G4CMPLewinSmithNIEL::PartitionNIEL(energy, material, Zin, Ain);
  }

  if (EmpEDepK && (EmpElow == EmpEhigh)) {
      
      G4Exception("G4CMPEmpiricalNIEL", "G4CMP1003", JustWarning,
          "Elow and Ehigh are equal. Defaulting to Lewin-Smith model.");
      useLewinSmith = true;
      
      return G4CMPLewinSmithNIEL::PartitionNIEL(energy, material, Zin, Ain);
  }

  G4double k = EmpkFixed;
  if (EmpEDepK) {
      G4double dk = Empkhigh - Empklow;
      G4double dE = EmpEhigh - EmpElow;
      k = Empklow + (energy-EmpElow)*dk/dE;
  }

  const G4double Z=GetEffectiveZ(material), A=GetEffectiveA(material)/(g/mole);
  G4Pow* g4pow = G4Pow::GetInstance();		// Fast, tabulated x^(2/3) etc.
  G4double z23 = g4pow->Z23(Z);

  G4double epsilon = 11.5/keV * energy / (Z*z23*z23);	// * Z^(-7/3)
  
  G4double h = (0.7*g4pow->powA(epsilon,0.6) + 3.*g4pow->powA(epsilon,0.15) + epsilon);

  return (k*h / (1.+k*h));
  }
}
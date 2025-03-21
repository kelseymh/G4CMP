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

// Constructor with debug prints
G4CMPEmpiricalNIEL::G4CMPEmpiricalNIEL() 
    : Empklow(DBL_MIN), Empkhigh(DBL_MIN), 
      EmpElow(DBL_MIN), EmpEhigh(DBL_MIN), 
      EmpkFixed(DBL_MIN), EmpEDepK(true) {

    G4CMPConfigManager* config = G4CMPConfigManager::Instance();

    // Fetch values from G4CMPConfigManager or set defaults
    Empklow = config->GetEmpklow();
    Empkhigh = config->GetEmpkhigh();
    EmpElow = config->GetEmpElow();
    EmpEhigh = config->GetEmpEhigh();
    EmpkFixed = config->GetEmpkFixed();
    EmpEDepK = config->GetEmpEDepK();
}

G4double G4CMPEmpiricalNIEL::
PartitionNIEL(G4double energy, const G4Material *material, G4double Zin, G4double Ain) const {
    if (!material) {
        G4Exception("G4CMPEmpiricalNIEL", "G4CMP1000", FatalErrorInArgument,
            "No material passed to partition function");
        return 1.;  // Won't get here after FATAL error
    }

    if (useLewinSmith) {
        G4cout << "Using Lewin-Smith model" << G4endl;
        return G4CMPLewinSmithNIEL::PartitionNIEL(energy, material, Zin, Ain);
    }

    if (EmpEDepK && (energy > EmpEhigh)) {
        if (firstCall_E) {
            G4Exception("G4CMPEmpiricalNIEL", "G4CMP1002", JustWarning,
                "Energy is greater than Ehigh. Lewin-Smith model is used.");
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
        k = Empklow + (energy - EmpElow) * dk / dE;
    }


    // Compute effective Z and A
    const G4double Z = GetEffectiveZ(material);
    const G4double A = GetEffectiveA(material) / (g/mole);

    // Compute epsilon and h
    G4Pow* g4pow = G4Pow::GetInstance();
    G4double z23 = g4pow->Z23(Z);
    G4double epsilon = 11.5 / keV * energy / (Z * z23 * z23);
    G4double h = (0.7 * g4pow->powA(epsilon, 0.6) + 3.0 * g4pow->powA(epsilon, 0.15) + epsilon);

    G4double result = (k * h / (1.0 + k * h));

    return result;
}

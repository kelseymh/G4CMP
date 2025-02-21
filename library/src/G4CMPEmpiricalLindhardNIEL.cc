/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPEmpiricalLindhardNIEL.cc
/// \brief NIEL calculation from Ionization yield measurement 
/// in a germanium CDMSlite detector using photo-neutron sources.
///
/// Computation of NIEL according to Lindhard (Lewin Smith) yield with 
/// a parameter k that depends on the energy.
///
/// Paper DOI: https://doi.org/10.1103/PhysRevD.105.122002
//
// 20250212  David Sadek

#include "globals.hh"
#include "G4CMPEmpiricalLindhardNIEL.hh"
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


G4CMPEmpiricalLindhardNIEL::G4CMPEmpiricalLindhardNIEL() 
    : Empiricalklow(DBL_MIN), Empiricalkhigh(DBL_MIN), 
      EmpiricalElow(DBL_MIN), EmpiricalEhigh(DBL_MIN), 
      EmpiricalkFixed(DBL_MIN), EmpiricalEnergyDependentK(true) {

    G4CMPConfigManager* config = G4CMPConfigManager::Instance();

    // Fetch values from G4CMPConfigManager or set defaults
    Empiricalklow = (config->GetEmpiricalklow() != DBL_MIN) ? config->GetEmpiricalklow() : 0.040;
    Empiricalkhigh = (config->GetEmpiricalkhigh() != DBL_MIN) ? config->GetEmpiricalkhigh() : 0.142;
    EmpiricalElow = (config->GetEmpiricalElow() != DBL_MIN) ? config->GetEmpiricalElow() : 0.39 * keV;
    EmpiricalEhigh = (config->GetEmpiricalEhigh() != DBL_MIN) ? config->GetEmpiricalEhigh() : 7.0 * keV;
    EmpiricalkFixed = (config->GetEmpiricalkFixed() != DBL_MIN) ? config->GetEmpiricalkFixed() : 0.158;
    EmpiricalEnergyDependentK = config->GetEmpiricalEnergyDependentK();
}

G4double G4CMPEmpiricalLindhardNIEL::
PartitionNIEL(G4double energy, const G4Material *material, G4double Zin,
		G4double Ain) const {

  if (!material) {
    G4Exception("G4CMPEmpiricalLindhardNIEL", "G4CMP1000", FatalErrorInArgument,
		  "No material passed to partition function");
    return 1.;		// Won't get here after FATAL error
  }
    
  if (firstCall) {
        	G4Exception("G4CMPEmpiricalLindhardNIEL", "G4CMP1001", JustWarning,
		    "Lindhard ionization yield model where the parameter k depends\n"
            "on the recoil energy is selected. Default settings are for Ge.\n"
            "If klow, khigh, Elow, Ehigh are not passed, default values from\n"
            "Phys. Rev. D 105, 122002 will be used.");
            firstCall = false;
  }
    
  if (EmpiricalEnergyDependentK && ((energy < EmpiricalElow) || (energy > EmpiricalEhigh))) {
      if (firstCall_E) {
          G4Exception("G4CMPEmpiricalLindhardNIEL", "G4CMP1002", JustWarning,
		    "Energy is out of bounds. Lewin-Smith model is used.");
         firstCall_E = false;
      }
      return G4CMPLewinSmithNIEL::PartitionNIEL(energy, material);
  }

  if (EmpiricalEnergyDependentK && (EmpiricalElow == EmpiricalEhigh))) {
      
      G4Exception("G4CMPEmpiricalLindhardNIEL", "G4CMP1003", JustWarning,
          "Elow and Ehigh are equal. Defaulting to Lewin-Smith model.");
      
      return G4CMPLewinSmithNIEL::PartitionNIEL(energy, material);
  }

  const G4double Z=GetEffectiveZ(material), A=GetEffectiveA(material)/(g/mole);
  G4Pow* g4pow = G4Pow::GetInstance();		// Fast, tabulated x^(2/3) etc.
  G4double z23 = g4pow->Z23(Z);

  G4double epsilon = 11.5/keV * energy / (Z*z23*z23);	// * Z^(-7/3)
  G4double k = EmpiricalEnergyDependentK ? Empiricalklow + ((Empiricalkhigh - Empiricalklow) / (EmpiricalEhigh - EmpiricalElow)) * (energy - EmpiricalElow) : EmpiricalkFixed;
  
  G4double h = (0.7*g4pow->powA(epsilon,0.6) + 3.*g4pow->powA(epsilon,0.15) + epsilon);

  return (k*h / (1.+k*h));
}


// Setter Functions
void G4CMPEmpiricalLindhardNIEL::SetEmpiricalklow(G4double kl) { Empiricalklow = kl; G4CMPConfigManager::Instance()->SetEmpiricalklow(kl); }
void G4CMPEmpiricalLindhardNIEL::SetEmpiricalkhigh(G4double kh) { Empiricalkhigh = kh; G4CMPConfigManager::Instance()->SetEmpiricalkhigh(kh); }
void G4CMPEmpiricalLindhardNIEL::SetEmpiricalElow(G4double el) { EmpiricalElow = el; G4CMPConfigManager::Instance()->SetEmpiricalElow(el); }
void G4CMPEmpiricalLindhardNIEL::SetEmpiricalEhigh(G4double eh) { EmpiricalEhigh = eh; G4CMPConfigManager::Instance()->SetEmpiricalEhigh(eh); }
void G4CMPEmpiricalLindhardNIEL::SetEmpiricalKFixed(G4double kf) { EmpiricalkFixed = kf; G4CMPConfigManager::Instance()->SetEmpiricalkFixed(kf); }
void G4CMPEmpiricalLindhardNIEL::SetEmpiricalEnergyDependentK(bool use) { EmpiricalEnergyDependentK = use; G4CMPConfigManager::Instance()->SetEmpiricalEnergyDependentK(use); }

// Getter Functions
G4double G4CMPEmpiricalLindhardNIEL::GetEmpiricalklow() const { return Empiricalklow; }
G4double G4CMPEmpiricalLindhardNIEL::GetEmpiricalkhigh() const { return Empiricalkhigh; }
G4double G4CMPEmpiricalLindhardNIEL::GetEmpiricalElow() const { return EmpiricalElow; }
G4double G4CMPEmpiricalLindhardNIEL::GetEmpiricalEhigh() const { return EmpiricalEhigh; }
G4double G4CMPEmpiricalLindhardNIEL::GetEmpiricalkFixed() const { return EmpiricalkFixed; }
bool G4CMPEmpiricalLindhardNIEL::GetEmpiricalEnergyDependentK() const { return EmpiricalEnergyDependentK; }
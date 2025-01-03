/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPSarkisNIEL.hh
/// \brief Non-ionizing energy loss calculation from Sarkis 2023.
///
/// Paper: https://journals.aps.org/pra/abstract/10.1103/PhysRevA.107.062811.
///
/// This ionization model was obtained from the Sarkis paper referenced above 
/// for Silicon ONLY so it doos not have (Z,A) dependence. The code will check
/// the effective Z of the input material. If the effZ is within
/// +/-1 of Silicon Z and A, the Sarkis model will be used, else,
/// Lindhard (LewinSmith) model will be used for NIEL calculations.
///
/// The model is obtained in the range of ~50 eV to 3 MeV above which
/// Lindahard (LewinSmith) model will be used for NIEL calculations

// 20230721  David Sadek - University of Florida (david.sadek@ufl.edu)
// 20240416  S. Zatschler -- Remove unused const A
// 20240705  M. Kelsey -- Restore A, and use it in test to check vs. SiA
// 20250102  M. Kelsey -- Fix warning messages to report correct name; move
//	       constructor here, to load data file once per job.

#include "globals.hh"
#include "G4CMPSarkisNIEL.hh"
#include "G4CMPConfigManager.hh"
#include "G4ExceptionSeverity.hh"
#include "G4SystemOfUnits.hh"
#include <fstream>
#include <algorithm>
#include <cmath>


// Constructor loads vector of data points for interpolation

G4CMPSarkisNIEL::G4CMPSarkisNIEL() {
  G4String fPath =
    G4CMPConfigManager::GetLatticeDir()+"/NIEL/SarkisSiIonYield.txt";
  std::ifstream inputFile(fPath);       // open the data file

  if (!inputFile.good()) {
    G4Exception("G4CMPSarkisNIEL", "G4CMP1101", FatalException,
		"Unable to open SarkisSiIonYield data file");
    return;
  }

  if (!lVector.Retrieve(inputFile, true)) {   // load the data from text file
    G4Exception("G4CMPSarkisNIEL", "G4CMP1102", FatalException,
		"Unable to successfully read SarkisSiIonYield data file");
    return;
  }

  lVector.SetSpline(true);		// Ensure that we can interpolate
  lVector.FillSecondDerivatives();

  G4cout << "G4CMPSarkisNIEL constructor:" << G4endl
	 << "lVector has " << lVector.GetVectorLength() << " entries,"
	 << " min energy " << lVector.GetLowEdgeEnergy(0)
	 << " max energy " << lVector.GetMaxEnergy() << " (no units!)"
	 << G4endl;
}

G4double G4CMPSarkisNIEL::
PartitionNIEL(G4double energy, const G4Material *material,G4double Zin,
		G4double Ain) const {
  if (!material) {
    G4Exception("G4CMPSarkisNIEL", "G4CMP1100", FatalErrorInArgument,
		  "No material passed to partition function");
    return 1.;		// Won't get here after FATAL error
  }

  // Get effective Z,A of the material
  const G4double Z = GetEffectiveZ(material);
  const G4double A = GetEffectiveA(material);
    
  // Check if the material is silicon or similar
  if (std::abs(Z-SiZ) < 0.5 && std::abs(A/(g/mole)-SiA) < 1.) {
    // Sarkis model below 3 MeV 
    if (energy <= 3*MeV) {
      return lVector.Value(energy, idx);	// Interpolate to get yield
    } else {
      // Lindhard model above 3 MeV
      if (firstCall) {
	G4Exception("G4CMPSarkisNIEL", "G4CMP1105", JustWarning,
		    "Sarkis model is obtained in the range of 50 eV to 3 MeV.\nAbove 3 MeV, Lindhard model will be used.");
	firstCall = false;
      }
      return G4CMPLewinSmithNIEL::PartitionNIEL(energy, material, Zin, Ain);
    }
  } else {
    if (firstCall) {
      G4Exception("G4CMPSarkisNIEL", "G4CMP1104", JustWarning,
		  "The input material is not Silicon.\nThe Lindhard model will be used for NIEL calculation.");
      firstCall = false;
    }
    return G4CMPLewinSmithNIEL::PartitionNIEL(energy, material, Zin, Ain);
  }
}

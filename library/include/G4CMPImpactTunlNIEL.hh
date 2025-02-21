/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPImpactTUNLNIEL.hh
/// \brief Non-ionizing energy loss calculation from IMPACT@TUNL 2023.
///
/// Computation of NIEL using the Emp model extracted from the
/// IMPACT@TUNL ionization yield measurements.  Link to the
/// paper:https://arxiv.org/abs/2303.02196.
///
/// This ionization model was obtained from the ionization yield
/// measurements in Silicon ONLY and it deos not have (Z,A)
/// dependence. The code will check the effective Z and A of the input
/// material, the effZ and effA are within +/-1 of Silicon Z and A, the
/// Impact model will be used, else, Lindhard(LewinSmith) model will be
/// used for NIEL calculations.
///
/// The model is obtained in the range of 100 eV to 10 keV. Above 10
/// keV, Lindhard(LewinSmith) model will be used.
//
// 20230721  David Sadek - University of Florida (david.sadek@ufl.edu)
// 20250102  M. Kelsey -- SiA should be in amu's (g/mol)

#ifndef G4CMPImpactTunlNIEL_hh
#define G4CMPImpactTunlNIEL_hh 1

#include "G4CMPLewinSmithNIEL.hh"


class G4CMPImpactTunlNIEL : public G4CMPLewinSmithNIEL {
public:
  G4CMPImpactTunlNIEL() {;}
  virtual ~G4CMPImpactTunlNIEL() {;}
  
  // return the fraction of the specified energy which will be deposited as NIEL
  // if an incoming particle with z1, a1 is stopped in the specified material
  // a1 is in atomic mass units, energy in native G4 energy units.
  //
  virtual G4double 
  PartitionNIEL(G4double energy, const G4Material *material, G4double Zin=0.,
		G4double Ain=0.) const override;

private:
  // A least-square fit is applied to the results on the ring detectors with 
  // an Emply chosen power-law function Y(Er)=Y10(Er/10000)^B

  const G4double B = 0.261;	// Best fit value +0.017, -0.011
  const G4double Y10 = 0.302;	// yield at 10 keV
  const G4double SiZ = 14.0;	// Z of Silicon
  const G4double SiA = 28.09;	// A of Silicon

  mutable bool firstCall = true; 	// Print warning messages only once
};

#endif	/* G4CMPImpactTunlNIEL_hh */

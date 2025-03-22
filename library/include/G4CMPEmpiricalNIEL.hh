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
// 20250212  David Sadek

#ifndef G4CMPEmpiricalNIEL_hh
#define G4CMPEmpiricalNIEL_hh 1

#include "G4VNIELPartition.hh"
#include "G4CMPLewinSmithNIEL.hh"
#include "G4SystemOfUnits.hh"
#include <cfloat>


class G4CMPEmpiricalNIEL : public G4CMPLewinSmithNIEL {
public:
    // Default constructor
    G4CMPEmpiricalNIEL();

    virtual ~G4CMPEmpiricalNIEL() override {;}

  // return the fraction of the specified energy which will be deposited as NIEL
  // if an incoming particle is stopped in the specified material,
  // energy in native G4 energy units.
  //
  // The parameter k is a function of energy

  // Setters
  void SetEmpklow(G4double kl) { Empklow = kl; }
  void SetEmpkhigh(G4double kh) { Empkhigh = kh; }
  void SetEmpElow(G4double el) { EmpElow = el; }
  void SetEmpEhigh(G4double eh) { EmpEhigh = eh; }
  void SetEmpkFixed(G4double fixedK) { EmpkFixed = fixedK; }
  void SetEmpEDepK(bool use) { EmpEDepK = use; }

  // Getters
  G4double GetEmpklow() const { return Empklow; }
  G4double GetEmpkhigh() const { return Empkhigh; }
  G4double GetEmpElow() const { return EmpElow; }
  G4double GetEmpEhigh() const { return EmpEhigh; }
  G4double GetEmpkFixed() const { return EmpkFixed; }
  bool GetEmpEDepK() const { return EmpEDepK; }

  virtual G4double 
  PartitionNIEL(G4double energy, const G4Material *material, G4double Zin=0.,
		G4double Ain=0.) const override;

private:

  // The following parameters are used to set the K vs E relation
    // Model fit parameters
  G4double Empklow;             // Set with /g4cmp/NIELPartition/Empirical/Empklow      value
  G4double Empkhigh;            // Set with /g4cmp/NIELPartition/Empirical/Empkhigh     value 
    // Model validity energy range
  G4double EmpElow;             // Set with /g4cmp/NIELPartition/Empirical/EmpElow      value units
  G4double EmpEhigh;            // Set with /g4cmp/NIELPartition/Empirical/EmpEhigh     value units
    // Flag to use Empirical Lindhard with energy-dependent k
  G4bool EmpEDepK;           // Set with /g4cmp/NIELPartition/Empirical/EmpEDepK  [true|false] 
    // If k is not energy dependent, provide/use kFixed
  G4double EmpkFixed;           // Set with /g4cmp/NIELPartition/Empirical/EmpkFixed    value

  mutable bool useLewinSmith = false; // if Elow=Ehigh, the Emp model blows up. 
  mutable bool firstCall_E = true; 	// Print warning messages only once
};

#endif	/* G4CMPEmpiricalNIEL_hh */

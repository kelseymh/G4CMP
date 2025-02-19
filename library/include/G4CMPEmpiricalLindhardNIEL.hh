/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPEmpiricalLindhardNIEL.hh
/// \brief NIEL calculation from Ionization yield measurement 
/// in a germanium CDMSlite detector using photo-neutron sources.
///
/// Computation of NIEL according to Lindhard (Lewin Smith) yield with 
/// a parameter k that depends on the energy.
///
/// Paper DOI: https://doi.org/10.1103/PhysRevD.105.122002
//
// 20250212  David Sadek

#ifndef G4CMPEmpiricalLindhardNIEL_hh
#define G4CMPEmpiricalLindhardNIEL_hh 1

#include "G4VNIELPartition.hh"
#include "G4CMPLewinSmithNIEL.hh"
#include "G4SystemOfUnits.hh"
#include <cfloat>


class G4CMPEmpiricalLindhardNIEL : public G4CMPLewinSmithNIEL {
public:
    // Default constructor
    G4CMPEmpiricalLindhardNIEL();

    virtual ~G4CMPEmpiricalLindhardNIEL() override {;}

  
  // return the fraction of the specified energy which will be deposited as NIEL
  // if an incoming particle is stopped in the specified material,
  // energy in native G4 energy units.
  //
  // The parameter k is a function of energy

  // Setters
  void SetEmpiricalklow(G4double kl);
  void SetEmpiricalkhigh(G4double kh);
  void SetEmpiricalElow(G4double el);
  void SetEmpiricalEhigh(G4double eh);
  void SetEmpiricalkFixed(G4double fixedK);
  void SetEmpiricalEnergyDependentK(bool use);


  // Getters
  G4double GetEmpiricalklow() const;
  G4double GetEmpiricalkhigh() const;
  G4double GetEmpiricalElow() const;
  G4double GetEmpiricalEhigh() const;
  G4double GetEmpiricalkFixed() const;
  bool GetEmpiricalEnergyDependentK() const;

  virtual G4double 
  PartitionNIEL(G4double energy, const G4Material *material, G4double Zin=0.,
		G4double Ain=0.) const override;

private:

  G4double Empiricalklow;      // Lower value of k
  G4double Empiricalkhigh;     // High value of k
  G4double EmpiricalElow;      // Lowest energy value 
  G4double EmpiricalEhigh;     // Highest energy value
  G4double EmpiricalkFixed;    // Value of k if not energy dependent nor passed
  bool EmpiricalEnergyDependentK;         // Flag to determine whether to use energy dependent k
  mutable bool firstCall = true; 	// Print warning messages only once
  mutable bool firstCall_E = true; 	// Print warning messages only once


};

#endif	/* G4CMPEmpiricalLindhardNIEL_hh */

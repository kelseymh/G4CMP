/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
// class G4CMPEqEMField,  created from G4EqMagElectricField by D. Brandt
//
// Class description:
//
// This is the right-hand side of equation of motion in a combined
// electric and magnetic field.
//
// 20140318  Need electron mass tensor as well as valley
// 20140331  Inherit from G4EqMagElectricField to handle holes as well as
//	     electrons.  Do local/global transformations; take valley index
//	     run-time configuration argument.
// 20140404  Drop unnecessary data members, using functions in G4LatticePhysical
// 20170525  Add default "rule of five" copy/move operators
// 20210920  Add verbosity with access to be used by G4CMPFieldManager

#ifndef G4CMPEqEMField_hh
#define G4CMPEqEMField_hh

#include "G4EqMagElectricField.hh"
#include "G4AffineTransform.hh"
#include "G4LatticePhysical.hh"
#include "G4RotationMatrix.hh"

class G4ElectroMagneticField;


class G4CMPEqEMField : public G4EqMagElectricField
{
public:
  G4CMPEqEMField(G4ElectroMagneticField *emField,
		 const G4LatticePhysical* lattice=0);

  virtual ~G4CMPEqEMField() {;} 

  G4CMPEqEMField(const G4CMPEqEMField&) = default;
  G4CMPEqEMField(G4CMPEqEMField&&) = default;
  G4CMPEqEMField& operator=(const G4CMPEqEMField&) = default;
  G4CMPEqEMField& operator=(G4CMPEqEMField&&) = default;

  // Turn on diagnostic messages (warning, this makes huge output files!)
  void SetVerboseLevel(G4int vb) { verboseLevel = vb; }

  // Replace physical lattice if track has changed volumes
  // NOTE:  Returns TRUE if lattice was actually changed
  G4bool ChangeLattice(const G4LatticePhysical* lattice);

  // Configure for local coordinates and electron valley axis
  void SetTransforms(const G4AffineTransform& lToG);
  void SetValley(size_t ivalley);
  void SetNoValley() { valleyIndex = -1; }		// Use this for holes

  // Configuration function from base class
  // NOTE: change of signature with G4 10.0
  virtual void SetChargeMomentumMass(G4ChargeState particleCharge,
				     G4double MomentumXc,
				     G4double mass);
  
  void EvaluateRhsGivenB(const G4double y[],
			 const G4double field[],
			 G4double dydx[]) const;
  // Given the value of the electromagnetic field, this function 
  // calculates the value of the derivative dydx.
  
private:
  const G4LatticePhysical* theLattice;
  G4int verboseLevel;			// For diagnostic messages

  G4double fCharge;	       		// Same as base class fElectrMagCof
  G4double fMass;	       		// Same as base class fElectrMagCof
  G4int valleyIndex;			// Index of current valley (-1 if none)

  G4AffineTransform fLocalToGlobal;	// Local vs. global coordinates
  G4AffineTransform fGlobalToLocal;
};

#endif

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

#ifndef G4CMPEqEMField_hh
#define G4CMPEqEMField_hh

#include "G4EqMagElectricField.hh"
#include "G4AffineTransform.hh"
#include "G4LatticePhysical.hh"
#include "G4RotationMatrix.hh"
#include "G4Version.hh"

class G4ElectroMagneticField;


class G4CMPEqEMField : public G4EqMagElectricField
{
public:
  G4CMPEqEMField(G4ElectroMagneticField *emField,
		 const G4LatticePhysical* lattice=0);

  ~G4CMPEqEMField() {;} 

  // Replace physical lattice if track has changed volumes
  // NOTE:  Returns TRUE if lattice was actually changed
  G4bool ChangeLattice(const G4LatticePhysical* lattice);

  // Configure for local coordinates and electron valley axis
  void SetTransforms(const G4AffineTransform& lToG);
  void SetValley(size_t ivalley);
  void SetNoValley();			// Use this for holes

  // Configuration function from base class
  virtual void SetChargeMomentumMass(G4double particleCharge, // in e+ units
				     G4double MomentumXc,
				     G4double mass);
  
#if G4VERSION_NUMBER >= 1000
  // Prevent the previous function from "hiding" this new base version
  virtual void SetChargeMomentumMass(G4ChargeState particleCharge,
				     G4double MomentumXc,
				     G4double mass) {
    SetChargeMomentumMass(particleCharge.GetCharge(), MomentumXc, mass);
  }
#endif
  
  void EvaluateRhsGivenB(const G4double y[],
			 const G4double field[],
			 G4double dydx[]) const;
  // Given the value of the electromagnetic field, this function 
  // calculates the value of the derivative dydx.
  
private:
  const G4LatticePhysical* theLattice;

  G4double fCharge;
  G4double fMass;

  G4AffineTransform fLocalToGlobal;	// Local vs. global coordinates
  G4AffineTransform fGlobalToLocal;

  void SetValleyTransform(const G4RotationMatrix& xform);
  G4bool useValley;			// Flag to avoid matrix op==()
  G4RotationMatrix normalToValley;	// Parameters for electron motion
  G4RotationMatrix valleyToNormal;
  G4RotationMatrix massInverse;
};

#endif

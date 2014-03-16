// $Id$
// class G4CMPEqEMField,  created from G4EqMagElectricField by D. Brandt
//
// Class description:
//
// This is the right-hand side of equation of motion in a combined
// electric and magnetic field.

#ifndef G4CMPEqEMField_hh
#define G4CMPEqEMField_hh

#include "G4EquationOfMotion.hh"
#include "G4AffineTransform.hh"
#include "G4Version.hh"

class G4ElectroMagneticField;


class G4CMPEqEMField : public G4EquationOfMotion
{
public:
  G4CMPEqEMField(G4ElectroMagneticField *emField,
		const G4AffineTransform& valleyXform);

  ~G4CMPEqEMField() {;} 

  void SetChargeMomentumMass(G4double particleCharge, // in e+ units
			     G4double MomentumXc,
			     G4double mass);
  
#if G4VERSION_NUMBER >= 1000
  // Prevent the previous function from "hiding" this new base version
  void SetChargeMomentumMass(G4ChargeState particleCharge,
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
  
  void SetValleyTransform(const G4AffineTransform& xform);
  inline G4AffineTransform GetValleyTransform() {return normalToValley;}
  
public:
  G4AffineTransform normalToValley;
  G4AffineTransform valleyToNormal;
  
private:
  G4double fElectroMagCof;
  G4double fMassCof;
};

#endif

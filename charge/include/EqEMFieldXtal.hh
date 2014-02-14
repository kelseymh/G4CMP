// CREATED FROM::
//$Id$
// GEANT4 tag $Name: geant4-09-03-patch-01 $
//
//
// class EqEMFieldXtal,  created from G4EqMagElectricField by D. Brandt
//
// Class description:
//
// This is the right-hand side of equation of motion in a combined
// electric and magnetic field.

#ifndef EQEMFIELDXTAL_hh
#define EQEMFIELDXTAL_hh

#include "G4EquationOfMotion.hh"
#include "G4AffineTransform.hh"
#include "G4Version.hh"

class G4ElectroMagneticField;


class EqEMFieldXtal : public G4EquationOfMotion
{
public:
  EqEMFieldXtal(G4ElectroMagneticField *emField );
  ~EqEMFieldXtal() {;} 

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

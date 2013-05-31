
// CREATED FROM::
//$Id: G4EqMagElectricField.hh,v 1.9 2006/06/29 18:22:03 gunter Exp $
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
#include "G4ElectroMagneticField.hh"

#include "G4AffineTransform.hh"
#include "G4ThreeVector.hh"

#define PI 3.1415926

class EqEMFieldXtal : public G4EquationOfMotion
{
  public:  // with description
  
  G4AffineTransform normalToValley;
  G4AffineTransform valleyToNormal;


    EqEMFieldXtal(G4ElectroMagneticField *emField )
      : G4EquationOfMotion( emField ) {
      normalToValley= G4AffineTransform(G4RotationMatrix(-PI/4, -PI/4,-PI/4));
      valleyToNormal= G4AffineTransform(G4RotationMatrix(-PI/4, -PI/4,-PI/4)).Inverse();
    }

    ~EqEMFieldXtal() {;} 

    void  SetChargeMomentumMass(G4double particleCharge, // in e+ units
                                G4double MomentumXc,
                                G4double mass);

    void EvaluateRhsGivenB(const G4double y[],
                           const G4double field[],
                                 G4double dydx[] ) const;
      // Given the value of the electromagnetic field, this function 
      // calculates the value of the derivative dydx.

  void SetValleyTransform(G4AffineTransform xform);

  private:

    G4double        fElectroMagCof ;
    G4double        fMassCof;
    G4ThreeVector 	me;
};

#endif

//
// MODIFIED::
// $Id$
// GEANT4 tag $Name: geant4-09-03-patch-01 $
//
// -------------------------------------------------------------------
// This version modified to allow oblique electron propagation
// D. Brandt
// 05 / 23 / 2011

#include "G4CMPEqEMField.hh"
#include "G4ElectroMagneticField.hh"
#include "G4PhysicalConstants.hh"


G4CMPEqEMField::G4CMPEqEMField(G4ElectroMagneticField *emField,
			       const G4RotationMatrix& valleyXform,
			       const G4RotationMatrix& electronMInv)
  : G4EquationOfMotion(emField) {
  SetValleyTransform(valleyXform);	// Default transform, may be replaced
  SetMassTensor(electronMInv);
}


void  
G4CMPEqEMField::SetChargeMomentumMass(G4double particleCharge, // e+ units
				     G4double /*particleMom*/,
				     G4double particleMass) {
   fElectroMagCof = electron_charge*particleCharge;
   fMassCof = particleMass*particleMass ; 
}


void G4CMPEqEMField::EvaluateRhsGivenB(const G4double y[],
				       const G4double field[],
				       G4double dydx[]) const { 
    G4ThreeVector pc(y[3], y[4], y[5]);
    G4ThreeVector p = pc/c_light;
    
    G4RotationMatrix mInv = valleyToNormal*massInv*normalToValley;

    G4ThreeVector v = mInv*p;
    G4double vel = v.mag();

    G4ThreeVector Efield(field[3], field[4], field[5]);
    G4ThreeVector retForce = electron_charge*c_light*Efield/vel;
    
    dydx[0] = v[0]/vel;		// Fill return buffer
    dydx[1] = v[1]/vel;
    dydx[2] = v[2]/vel;
    dydx[3] = retForce[0];
    dydx[4] = retForce[1];
    dydx[5] = retForce[2];
    dydx[6] = 0.;		//not used
    dydx[7] = 1./vel;		// Lab Time of flight
}


void G4CMPEqEMField::SetValleyTransform(const G4RotationMatrix& xform) {
    normalToValley = xform;
    valleyToNormal = xform.inverse();
}

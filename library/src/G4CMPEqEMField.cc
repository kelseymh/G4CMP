//
// MODIFIED::
// $Id$
// GEANT4 tag $Name: geant4-09-03-patch-01 $
//
//  10.11.98   V.Grichine
//
// -------------------------------------------------------------------
// This version modified to allow oblique electron propagation
// D. Brandt
// 05 / 23 / 2011

#include "G4CMPEqEMField.hh"
#include "G4ElectroMagneticField.hh"
#include "G4PhysicalConstants.hh"
#include "G4RotationMatrix.hh"


G4CMPEqEMField::G4CMPEqEMField(G4ElectroMagneticField *emField,
			     const G4AffineTransform& valleyXform)
  : G4EquationOfMotion(emField) {
  SetValleyTransform(valleyXform);	// Default transform, may be replaced
}


void  
G4CMPEqEMField::SetChargeMomentumMass(G4double particleCharge, // e+ units
				     G4double /*particleMom*/,
				     G4double particleMass)
{
   fElectroMagCof = electron_charge*particleCharge;
   fMassCof = particleMass*particleMass ; 
   /*
   G4cout<<"\nEqEMFIeldXtal::SetChargeMomentumMass: using normalToValley xform: \n\t"<<normalToValley.NetRotation().colX() <<"\n\t"<<normalToValley.NetRotation().colY()<<"\n\t"<<normalToValley.NetRotation().colZ();

   G4cout<<"\nEqEMFIeldXtal::SetChargeMomentumMass: using valleyToNormal xform: \n\t"<<valleyToNormal.NetRotation().colX() <<"\n\t"<<valleyToNormal.NetRotation().colY()<<"\n\t"<<valleyToNormal.NetRotation().colZ();

   G4cout<<"EqEMField::SetChargeMomentumMass: Using the affine xform, a "<<G4ThreeVector(1*mm,0,0)/mm<< " vector becomes: "<<valleyToNormal.TransformPoint(normalToValley.TransformPoint(G4ThreeVector(1*mm,0,0)))/mm;
   */
}


void
G4CMPEqEMField::EvaluateRhsGivenB(const G4double y[],
			                const G4double field[],
				              G4double dydx[] ) const
{ 
    G4ThreeVector pc = G4ThreeVector(y[3], y[4], y[5]);
    G4ThreeVector p = pc/c_light;
    
    G4double me = electron_mass_c2/c_squared;
    G4RotationMatrix mInv = 
      valleyToNormal.NetRotation()*G4Rep3x3(1/1.588/me,   0.0    , 0.0,
					    0.0     , 1/.081/me, 0.0, 
					    0.0     ,   0.0    , 1/.081/me)
      *normalToValley.NetRotation();

    G4ThreeVector v = mInv*p;
    
    dydx[0] = v[0]/v.mag();
    dydx[1] = v[1]/v.mag();
    dydx[2] = v[2]/v.mag();
    
    G4ThreeVector Efield = G4ThreeVector(field[3], field[4], field[5]);
    G4ThreeVector retForce = electron_charge*c_light*Efield/v.mag();
    
    dydx[3] = retForce[0];
    dydx[4] = retForce[1];
    dydx[5] = retForce[2];

    dydx[6] = 0.;//not used

    // Lab Time of flight
    G4double inverse_velocity = 1/v.mag();
    dydx[7] = inverse_velocity;
    return ;
}


void G4CMPEqEMField::SetValleyTransform(const G4AffineTransform& xform) {
    normalToValley = xform;
    valleyToNormal = normalToValley.Inverse();
}

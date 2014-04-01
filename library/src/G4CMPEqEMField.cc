//
// MODIFIED::
// $Id$
// GEANT4 tag $Name: geant4-09-03-patch-01 $
//
// -------------------------------------------------------------------
// This version modified to allow oblique electron propagation
// D. Brandt
// 05 / 23 / 2011
//
// 20140331  Inherit from G4EqMagElectricField to handle holes as well as
//	     electrons.  Do local/global transformations; take valley index
//	     run-time configuration argument.

#include "G4CMPEqEMField.hh"
#include "G4ElectroMagneticField.hh"
#include "G4PhysicalConstants.hh"


G4CMPEqEMField::G4CMPEqEMField(G4ElectroMagneticField *emField,
			       const G4LatticePhysical* lattice)
  : G4EqMagElectricField(emField), theLattice(lattice), fCharge(0.),
    useValley(false) {
  if (lattice) massInverse = lattice->GetMInvTensor();
}


// Replace physical lattice if track has changed volumes
// NOTE:  Returns TRUE if lattice was actually changed

G4bool G4CMPEqEMField::ChangeLattice(const G4LatticePhysical* lattice) {
  G4bool newLat = (lattice != theLattice);
  if (newLat) {
    theLattice = lattice;
    if (lattice) massInverse = lattice->GetMInvTensor();
  }

  return newLat;
}


// Specify local-to-global transformation before field call
// Typically, this will be called from FieldManager::ConfigureForTrack()

void G4CMPEqEMField::SetTransforms(const G4AffineTransform& lToG) {
  fGlobalToLocal = fLocalToGlobal = lToG;
  fGlobalToLocal.Invert();
}


// Specify which valley to use for electrons, or no valley at all

void G4CMPEqEMField::SetValley(size_t ivalley) {
  if (theLattice && ivalley<theLattice->NumberOfValleys()) {
    SetValleyTransform(theLattice->GetValley(ivalley));
    useValley = true;
  } else {
    SetNoValley();
  }
}

void G4CMPEqEMField::SetNoValley() {
  SetValleyTransform(G4RotationMatrix::IDENTITY);
  useValley = false;
}


void G4CMPEqEMField::SetValleyTransform(const G4RotationMatrix& xform) {
  normalToValley = xform;
  valleyToNormal = xform.inverse();
}


// Configuration function must call through to base class
// NOTE: change of signature with G4 10.0

#if G4VERSION_NUMBER >= 1000
void G4CMPEqEMField::SetChargeMomentumMass(G4ChargeState particleCharge,
					   G4double MomentumXc,
					   G4double mass) {
  G4EqMagElectricField::SetChargeMomentumMass(particleCharge, MomentumXc, mass);
  fCharge = particleCharge.GetCharge() * electron_charge;
}
#else
void G4CMPEqEMField::SetChargeMomentumMass(G4double particleCharge,
					   G4double MomentumXc,
					   G4double mass) {
  G4EqMagElectricField::SetChargeMomentumMass(particleCharge, MomentumXc, mass);
  fCharge = particleCharge * electron_charge;
}
#endif
  

// Field evaluation:  Given momentum (y) and field, return velocity, force

void G4CMPEqEMField::EvaluateRhsGivenB(const G4double y[],
				       const G4double field[],
				       G4double dydx[]) const {
  // No lattice behaviour, just use base class
  if (!useValley) {
    G4EqMagElectricField::EvaluateRhsGivenB(y, field, dydx);
#ifdef G4CMP_DEBUG
    G4cout << "Base class returns dydx { "
	   << dydx[0] << " , " << dydx[1] << " , " << dydx[2] << " , "
	   << dydx[3] << " , " << dydx[4] << " , " << dydx[5] << " , "
	   << dydx[6] << " , " << dydx[7] << " }" << G4endl;
#endif
    return;
  }

  G4ThreeVector pc(y[3], y[4], y[5]);	// Momentum in G4 "natural units"
  G4ThreeVector p = pc/c_light;		// Convert e.g., MeV to MeV/c

  // Momentum to velocity conversion must be done in local coordinates
  fGlobalToLocal.ApplyAxisTransform(p);

  G4RotationMatrix mInv = valleyToNormal*massInverse*normalToValley;
  G4ThreeVector v = mInv*p;
  fLocalToGlobal.ApplyAxisTransform(v);

  G4double vel = v.mag();
  
  G4ThreeVector Efield(field[3], field[4], field[5]);
  G4ThreeVector retForce = fCharge * c_light * Efield/vel;
  
  dydx[0] = v[0]/vel;		// Velocity direction
  dydx[1] = v[1]/vel;
  dydx[2] = v[2]/vel;
  dydx[3] = retForce[0];	// Applied force
  dydx[4] = retForce[1];
  dydx[5] = retForce[2];
  dydx[6] = 0.;			// not used
  dydx[7] = 1./vel;		// Lab Time of flight

#ifdef G4CMP_DEBUG
    G4cout << "Valley calc returns dydx { "
	   << dydx[0] << " , " << dydx[1] << " , " << dydx[2] << " , "
	   << dydx[3] << " , " << dydx[4] << " , " << dydx[5] << " , "
	   << dydx[6] << " , " << dydx[7] << " }" << G4endl;
#endif
 }

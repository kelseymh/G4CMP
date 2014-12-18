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
// 20140404  Drop unnecessary data members, using functions in G4LatticePhysical
// 20140501  Fix sign flip in electron charge calculation.
// 20141217  Avoid floating-point division by using vinv = 1/v.mag()

#include "G4CMPEqEMField.hh"
#include "G4ElectroMagneticField.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"


G4CMPEqEMField::G4CMPEqEMField(G4ElectroMagneticField *emField,
			       const G4LatticePhysical* lattice)
  : G4EqMagElectricField(emField), theLattice(lattice), fCharge(0.),
    valleyIndex(-1) {;}


// Replace physical lattice if track has changed volumes
// NOTE:  Returns TRUE if lattice was actually changed

G4bool G4CMPEqEMField::ChangeLattice(const G4LatticePhysical* lattice) {
  G4bool newLat = (lattice != theLattice);
  theLattice = lattice;
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
    valleyIndex = ivalley;
  } else {
    valleyIndex = -1;
  }
}


// Configuration function must call through to base class
// NOTE: change of signature with G4 10.0

#if G4VERSION_NUMBER >= 1000
void G4CMPEqEMField::SetChargeMomentumMass(G4ChargeState particleCharge,
					   G4double MomentumXc,
					   G4double mass) {
  G4EqMagElectricField::SetChargeMomentumMass(particleCharge, MomentumXc, mass);
  fCharge = particleCharge.GetCharge() * eplus;
}
#else
void G4CMPEqEMField::SetChargeMomentumMass(G4double particleCharge,
					   G4double MomentumXc,
					   G4double mass) {
  G4EqMagElectricField::SetChargeMomentumMass(particleCharge, MomentumXc, mass);
  fCharge = particleCharge * eplus;
}
#endif
  

// Field evaluation:  Given momentum (y) and field, return velocity, force

void G4CMPEqEMField::EvaluateRhsGivenB(const G4double y[],
				       const G4double field[],
				       G4double dydx[]) const {
  // No lattice behaviour, just use base class
  if (valleyIndex < 0) {
    G4EqMagElectricField::EvaluateRhsGivenB(y, field, dydx);
    return;
  }

  // Momentum to velocity conversion must be done in local coordinates
  G4ThreeVector p(y[3], y[4], y[5]);
  fGlobalToLocal.ApplyAxisTransform(p);

  G4ThreeVector v = theLattice->MapPtoV_el(valleyIndex, p);
  fLocalToGlobal.ApplyAxisTransform(v);
  G4double vinv = 1./v.mag();
  
  G4ThreeVector Efield(field[3], field[4], field[5]);
  G4ThreeVector retForce = fCharge * Efield * c_light*vinv;
  
  dydx[0] = v.x()*vinv;		// Velocity direction
  dydx[1] = v.y()*vinv;
  dydx[2] = v.z()*vinv;
  dydx[3] = retForce.x();	// Applied force
  dydx[4] = retForce.y();
  dydx[5] = retForce.z();
  dydx[6] = 0.;			// not used
  dydx[7] = vinv;		// Lab Time of flight (sec/mm)
}

/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
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
// 20150528  Add debugging output
// 20190802  Check if field is aligned or anti-aligned with valley, apply
//	     transform to valley axis "closest" to field direction.
// 20210921  Add detailed debugging output, protected with G4CMP_DEBUG flag

#include "G4CMPEqEMField.hh"
#include "G4CMPConfigManager.hh"
#include "G4ElectroMagneticField.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"


G4CMPEqEMField::G4CMPEqEMField(G4ElectroMagneticField *emField,
			       const G4LatticePhysical* lattice)
  : G4EqMagElectricField(emField), theLattice(lattice), 
    verboseLevel(G4CMPConfigManager::GetVerboseLevel()),
    fCharge(0.), fMass(0.), valleyIndex(-1) {;}


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

void G4CMPEqEMField::SetChargeMomentumMass(G4ChargeState particleCharge,
					   G4double MomentumXc,
					   G4double mass) {
  G4EqMagElectricField::SetChargeMomentumMass(particleCharge, MomentumXc, mass);
  fCharge = particleCharge.GetCharge() * eplus;
  fMass = mass/c_squared;
}


// Field evaluation:  Given momentum (y) and field, return velocity, force

void G4CMPEqEMField::EvaluateRhsGivenB(const G4double y[],
				       const G4double field[],
				       G4double dydx[]) const {
  // No lattice behaviour, just use base class
  if (valleyIndex == -1) {
    G4EqMagElectricField::EvaluateRhsGivenB(y, field, dydx);
    return;
  }

#ifdef G4CMP_DEBUG
  if (verboseLevel > 2) {
    G4cout << "G4CMPEqEMField"
	   << " @ (" << y[0] << "," << y[1] << "," << y[2] << ") mm" << G4endl
	   << " (q,m) " << fCharge/eplus << " e+ "
	   << fMass*c_squared/electron_mass_c2 << " m_e"
	   << " valley " << valleyIndex << G4endl;
  }
#endif

  /* This part is confusing. "Momentum" reported by G4 is not really the
   * momentum for charge carriers with valleys. It's just the velocity times
   * the defined scalar mass.
   *
   * So we need to calculate the true dp/dx, and then transform it into dv/dx
   * and then multiply that by the mass to get this "pseudomomentum."
   */
  G4ThreeVector v = G4ThreeVector(y[3], y[4], y[5])/fMass/c_light;
  G4double vinv = 1./v.mag();

#ifdef G4CMPDEBUG
  if (verboseLevel>2) {
    G4cout << " pc (" << y[3] << "," << y[4] << "," << y[5] << ") MeV" << G4endl
	   << " v " << v/(km/s) << " km/s " << " TOF " << vinv/(ns/mm)
	   << " ns/mm" << G4endl;
  }
#endif

  G4ThreeVector Efield(field[3], field[4], field[5]);
  G4ThreeVector force = fCharge * Efield;

#ifdef G4CMPDEBUG
  if (verboseLevel > 2) {
    G4cout << " E " << Efield/(volt/cm) << " " << Efield.mag()/(volt/cm)
	   << " V/cm" << G4endl
	   << " q*E " << force/(eV/m) << " eV/m" << G4endl;
  }
#endif

  fGlobalToLocal.ApplyAxisTransform(force);
  theLattice->RotateToLattice(force);

  // Since F is proportional to dp, it must transform like momentum.
  const G4RotationMatrix& vToN = theLattice->GetValley(valleyIndex);
  const G4RotationMatrix& nToV = theLattice->GetValleyInv(valleyIndex);

#ifdef G4CMPDEBUG
  if (verboseLevel > 2) {
    G4cout << " q*E (lattice) " << force/(eV/m) << G4endl
	   << " q*E (valley) " << vToN*force/(eV/m) << G4endl
	   << " q*E/m-tensor " << theLattice->GetMInvTensor()*(vToN*force)/(eV/m)
	   << G4endl;
  }
#endif

  force = nToV*(theLattice->GetMInvTensor()*(vToN*force));
#ifdef G4CMPDEBUG
  if (verboseLevel > 2) G4cout << " q*E/m (lattice) " << force/(eV/m) << G4endl;
#endif

  force *= fMass * vinv * c_light;
  theLattice->RotateToSolid(force);

#ifdef G4CMPDEBUG
  if (verboseLevel > 2) {
    G4cout << " force (local) " << force/(eV/m) << " " << force.mag()/(eV/m)
	   << " eV/m" << G4endl;
  }
#endif

  // Restore effective force to global coordinates for G4Transporation
  fLocalToGlobal.ApplyAxisTransform(force);

#ifdef G4CMPDEBUG
  if (verboseLevel > 2) {
    G4cout << " force (global) " << force/(eV/m) << " eV/m" << G4endl;
  }
#endif

  dydx[0] = v.x()*vinv;		// Velocity direction
  dydx[1] = v.y()*vinv;
  dydx[2] = v.z()*vinv;
  dydx[3] = force.x();		// Applied force in H-V space
  dydx[4] = force.y();
  dydx[5] = force.z();
  dydx[6] = 0.;			// not used
  dydx[7] = vinv;		// Lab Time of flight (sec/mm)
}

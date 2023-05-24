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
// 20210922  Field transformation should be Herring-Vogt, with SqrtInvTensor.
// 20211004  Compute velocity from true momentum, local-to-global as needed.
//		Clarify field transform and force calculation.
// 20211007  Insert debugging output for each step of E-field transformation.
// 20211012  Apply scale factor to conserve energy averaged over many electrons

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
#ifdef G4CMP_DEBUG
  if (verboseLevel>2) {
    G4cout << "G4CMPEqEMField::SetChargeMomentumMass "
	   << particleCharge.GetCharge() << " " << MomentumXc << " MeV "
	   << mass/electron_mass_c2 << " m_e" << G4endl;
  }
#endif

  G4EqMagElectricField::SetChargeMomentumMass(particleCharge, MomentumXc, mass);
  fCharge = particleCharge.GetCharge() * eplus;
  fMass = mass;
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

  // Get kinematics into more usable form
  pos.set(y[0],y[1],y[2]);			// Position
  mom.set(y[3],y[4],y[5]);			// Momentum
  Efield.set(field[3],field[4],field[5]);	// Electric field
  G4double Emag = Efield.mag();

  force = Efield;	// Apply transforms here so Efield stays original

#ifdef G4CMP_DEBUG
  if (verboseLevel>2) {
    G4cout << "G4CMPEqEMField" << " @ " << pos << " mm" << G4endl
	   << " (q,m) " << fCharge/eplus << " e+ "
	   << fMass/electron_mass_c2 << " m_e"
	   << " valley " << valleyIndex << G4endl
	   << " pc " << mom << " " << mom.mag() << " MeV" << G4endl;
  }
#endif
  // G4double Energy = std::sqrt( y[3]*y[3]+y[4]*y[4]+y[5]*y[5] + fMass*fMass );
  // G4double MEff = theLattice->GetElectronEffectiveMass(valleyIndex, mom);
  // mom = mom.unit() * std::sqrt(Energy*(Energy+2.*fMass));
  /* "Momentum" reported by G4 is the true momentum.
   */

  vel = mom;
  vel *= c_light/(electron_mass_c2);		// v = p/m = pc/c / mc^2/c^2 = pc/(mc^2/c)
  G4double vinv = 1./vel.mag();

  momdir = vel.unit();

#ifdef G4CMP_DEBUG
  if (verboseLevel>2) {
    G4cout << " v " << vel/(km/s) << " " << vel.mag()/(km/s) << " km/s"
	   << G4endl << " TOF (1/v) " << vinv/(ns/mm) << " ns/mm"
	   << " c/v " << vinv*c_light << G4endl
	   << " E-field         " << Efield/(volt/cm) << " "
	   << Emag/(volt/cm) << " V/cm" << G4endl;
  }
#endif

  fGlobalToLocal.ApplyAxisTransform(force);
#ifdef G4CMP_DEBUG
  if (verboseLevel>2)
    G4cout << " Field (loc)     " << force/(volt/cm) << " "
	   << force.mag()/(volt/cm) << G4endl;
#endif

  theLattice->RotateToLattice(force);
#ifdef G4CMP_DEBUG
  if (verboseLevel>2)
    G4cout << " Field (lat)     " << force/(volt/cm) << " "
	   << force.mag()/(volt/cm) << G4endl;
#endif

  // Rotate force into and out of valley frame, applying Herring-Vogt transform
  const G4RotationMatrix& nToV = theLattice->GetValley(valleyIndex);
  const G4RotationMatrix& vToN = theLattice->GetValleyInv(valleyIndex);

  force.transform(nToV);			// Rotate to valley
#ifdef G4CMP_DEBUG
  if (verboseLevel>2)
    G4cout << " Field (val)     " << force/(volt/cm) << " "
	   << force.mag()/(volt/cm) << G4endl;
#endif

  force *= theLattice->GetMInvTensor();
  force *= electron_mass_c2/c_squared;
#ifdef G4CMP_DEBUG
  if (verboseLevel>2)
    G4cout << " Field (H-V)     " << force/(volt/cm) << " "
	   << force.mag()/(volt/cm) << G4endl;
#endif

  force.transform(vToN);			// Back to lattice
#ifdef G4CMP_DEBUG
  if (verboseLevel>2)
    G4cout << " Field (H-V lat) " << force/(volt/cm) << " "
	   << force.mag()/(volt/cm) << G4endl;
#endif

  theLattice->RotateToSolid(force);		// Back to crystal frame
#ifdef G4CMP_DEBUG
  if (verboseLevel>2)
    G4cout << " Field (H-V loc) " << force/(volt/cm) << " "
	   << force.mag()/(volt/cm) << G4endl;
#endif

  // Restore field to global coordinate frame for G4Transporation
  fLocalToGlobal.ApplyAxisTransform(force);
#ifdef G4CMP_DEBUG
  if (verboseLevel>2)
    G4cout << " Field (H-V glb) " << force/(volt/cm) << " "
	   << force.mag()/(volt/cm) << G4endl;
#endif

  // dp/ds = qE/beta
  force *= fCharge*vinv*c_light;;

#ifdef G4CMP_DEBUG
  if (verboseLevel>2) {
    G4cout << " q*Ec/v (scaled) " << force/(eV/m) << " " << force.mag()/(eV/m)
	   << " eV/m" << G4endl;
  }
#endif

  // Populate output buffer
  dydx[0] = momdir.x();		// Momentum direction
  dydx[1] = momdir.y();
  dydx[2] = momdir.z();
  dydx[3] = force.x();		// Effective force in H-V, global coordinates
  dydx[4] = force.y();
  dydx[5] = force.z();
  dydx[6] = 0.;			// not used
  dydx[7] = vinv;		// Lab Time of flight (ns/mm)
}

/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file materials/src/G4LatticePhysical.cc
/// \brief Implementation of the G4LatticePhysical class
//
// $Id$
//
// 20131115  Save rotation results in local variable, report verbosely
// 20131116  Replace G4Transform3D with G4RotationMatrix
// 20140319  Add output functions for diagnostics
// 20140321  Move placement transformations to G4CMPProcessUtils, put
//		lattice orientation into ctor arguments
// 20140401  Add valley momentum calculations
// 20140408  Move vally momentum calcs to G4LatticeLogical
// 20140425  Add "effective mass" calculation for electrons
// 20150601  Add mapping from electron velocity back to momentum
// 20160517  Replace unit vectors with CLHEP built-in values
// 20160608  Drop (theta,phi) lattice orientation function.
// 20170525  Drop empty destructor to allow default "rule of five" semantics
// 20170928  Replce "polarizationState" with "mode"
// 20190801  M. Kelsey -- Use G4ThreeVector buffer instead of pass-by-value
// 20200520  For MT thread safety, wrap G4ThreeVector buffer in function to
//		return thread-local instance.
// 20211021  Wrap verbose output in #ifdef G4CMP_DEBUG for performace
// 20220921  G4CMP-319 -- Add utilities for thermal (Maxwellian) distributions
// 20250507  G4CMP-480 -- Swap rotation matrix for local<-->lattice transforms.

#include "G4LatticePhysical.hh"
#include "G4CMPConfigManager.hh"
#include "G4LatticeLogical.hh"
#include "G4PhysicalConstants.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"


// Null vector defined for convenience (avoid memory churn)

namespace {
  G4ThreeVector nullVec(0,0,0);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Default constructor

G4LatticePhysical::G4LatticePhysical()
  : verboseLevel(0), fLattice(0), hMiller(0), kMiller(0), lMiller(0),
    fRot(0.), fTemperature(-1.) {;}

// Set lattice orientation (relative to G4VSolid) with Miller indices

G4LatticePhysical::G4LatticePhysical(const G4LatticeLogical* Lat,
				     G4int h, G4int k, G4int l, G4double rot)
  : verboseLevel(0), fLattice(Lat), fTemperature(-1.) {
  SetMillerOrientation(h, k, l, rot);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Align Miller normal vector (hkl) with +Z axis, and rotation about axis
void G4LatticePhysical::SetMillerOrientation(G4int h, G4int k, G4int l,
					     G4double rot) {
  if (verboseLevel) {
    G4cout << "G4LatticePhysical::SetMillerOrientation(" << h << " "
	   << k << " " << l << ", " << rot/deg << " deg)" << G4endl;
  }

  hMiller = h;
  kMiller = k;
  lMiller = l;
  fRot = rot;

  G4ThreeVector norm = (h*GetBasis(0)+k*GetBasis(1)+l*GetBasis(2)).unit();

  if (verboseLevel>1) G4cout << " norm = " << norm << G4endl;

  // Aligns geometry +Z axis with lattice (hkl) normal
  fOrient = G4RotationMatrix::IDENTITY;
  fOrient.rotateZ(rot).rotateY(norm.theta()).rotateZ(norm.phi());
  fInverse = fOrient.inverse();

  if (verboseLevel>1) G4cout << " fOrient = " << fOrient << G4endl;

  // FIXME:  Is this equivalent to (phi,theta,rot) Euler angles???
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Return temperature assigned to lattice/volume, or global parameter

G4double G4LatticePhysical::GetTemperature() const {
  return (fTemperature < 0. ? G4CMPConfigManager::GetTemperature()
	  : fTemperature);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Rotate input vector between lattice and solid orientations

const G4ThreeVector&
G4LatticePhysical::RotateToLattice(G4ThreeVector& dir) const {
  return dir.transform(fInverse);
}

const G4ThreeVector& 
G4LatticePhysical::RotateToSolid(G4ThreeVector& dir) const {
  return dir.transform(fOrient);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

///////////////////////////////
//Loads the group velocity in m/s
/////////////////////////////
G4double G4LatticePhysical::MapKtoV(G4int mode, const G4ThreeVector& k) const {
#ifdef G4CMP_DEBUG
  if (verboseLevel>1) G4cout << "G4LatticePhysical::MapKtoV " << k << G4endl;
#endif

  RotateToLattice(tempvec()=k);
#ifdef G4CMP_DEBUG
  if (verboseLevel>1) G4cout << " in lattice frame " << tempvec() << G4endl;
#endif

  return fLattice->MapKtoV(mode, tempvec());
}

///////////////////////////////
//Loads the normalized direction vector along VG
///////////////////////////////
G4ThreeVector G4LatticePhysical::MapKtoVDir(G4int mode, const G4ThreeVector& k) const {
#ifdef G4CMP_DEBUG
  if (verboseLevel>1) G4cout << "G4LatticePhysical::MapKtoVDir " << k << G4endl;
#endif

  RotateToLattice(tempvec()=k);
#ifdef G4CMP_DEBUG
  if (verboseLevel>1) G4cout << " in lattice frame " << tempvec() << G4endl;
#endif

  G4ThreeVector VG = fLattice->MapKtoVDir(mode, tempvec());
#ifdef G4CMP_DEBUG
  if (verboseLevel>1) G4cout << " VDir (lattice) " << tempvec() << G4endl;
#endif

  return RotateToSolid(VG);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector
G4LatticePhysical::MapEkintoP(G4int iv, const G4ThreeVector& pdir, const G4double Ekin) const {
#ifdef G4CMP_DEBUG
  if (verboseLevel>1)
    G4cout << "G4LatticePhysical::MapEkintoP " << iv << " " << pdir << " " << Ekin << G4endl;
#endif

  RotateToLattice(tempvec()=pdir);
#ifdef G4CMP_DEBUG
  if (verboseLevel>1) G4cout << " in lattice frame " << tempvec() << G4endl;
#endif

  G4ThreeVector p = fLattice->MapEkintoP(iv, tempvec(), Ekin);
#ifdef G4CMP_DEBUG
  if (verboseLevel>1) G4cout << " P (lattice) " << tempvec() << G4endl;
#endif

  return RotateToSolid(p);
}

G4double G4LatticePhysical::MapPtoEkin(G4int iv, const G4ThreeVector& p) const {
#ifdef G4CMP_DEBUG
  if (verboseLevel>1)
    G4cout << "G4LatticePhysical::MapPtoEkin " << iv << " " << p << G4endl;
#endif

  RotateToLattice(tempvec()=p);
#ifdef G4CMP_DEBUG
  if (verboseLevel>1) {
    G4cout << " in lattice frame " << tempvec() << G4endl
	   << " returning Ekin " << fLattice->MapPtoEkin(iv, tempvec())
	   << G4endl;
  }
#endif

  return fLattice->MapPtoEkin(iv, tempvec());
}

G4double G4LatticePhysical::MapV_elToEkin(G4int iv, const G4ThreeVector& v) const {
#ifdef G4CMP_DEBUG
  if (verboseLevel>1)
    G4cout << "G4LatticePhysical::MapV_elToEkin " << iv << " " << v << G4endl;
#endif

  RotateToLattice(tempvec()=v);
#ifdef G4CMP_DEBUG
  if (verboseLevel>1) {
    G4cout << " in lattice frame " << tempvec() << G4endl
	   << " returning Ekin " << fLattice->MapV_elToEkin(iv, tempvec())
	   << G4endl;
  }
#endif

  return fLattice->MapV_elToEkin(iv, tempvec());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Convert electron momentum to valley velocity, wavevector, and HV vector

G4ThreeVector 
G4LatticePhysical::MapPtoV_el(G4int ivalley, const G4ThreeVector& p_e) const {
#ifdef G4CMP_DEBUG
  if (verboseLevel>1)
    G4cout << "G4LatticePhysical::MapPtoV_el " << ivalley << " " << p_e
	   << G4endl;
#endif

  RotateToLattice(tempvec()=p_e);
#ifdef G4CMP_DEBUG
  if (verboseLevel>1) G4cout << " in lattice frame " << tempvec() << G4endl;
#endif

  tempvec() = fLattice->MapPtoV_el(ivalley, tempvec());
#ifdef G4CMP_DEBUG
  if (verboseLevel>1) G4cout << " V_el (lattice) " << tempvec() << G4endl;
#endif

  return RotateToSolid(tempvec());
}

G4ThreeVector 
G4LatticePhysical::MapV_elToP(G4int ivalley, const G4ThreeVector& v_e) const {
#ifdef G4CMP_DEBUG
  if (verboseLevel>1)
    G4cout << "G4LatticePhysical::MapV_elToP " << ivalley << " " << v_e
	   << G4endl;
#endif

  RotateToLattice(tempvec()=v_e);
#ifdef G4CMP_DEBUG
  if (verboseLevel>1) G4cout << " in lattice frame " << tempvec() << G4endl;
#endif

  tempvec() = fLattice->MapV_elToP(ivalley, tempvec());
#ifdef G4CMP_DEBUG
  if (verboseLevel>1) G4cout << " p (lattice) " << tempvec() << G4endl;
#endif

  return RotateToSolid(tempvec());
}

G4ThreeVector 
G4LatticePhysical::MapPToP_Q(G4int ivalley, const G4ThreeVector& P) const {
#ifdef G4CMP_DEBUG
  if (verboseLevel>1)
    G4cout << "G4LatticePhysical::MapPToP_Q " << ivalley << " " << P
	   << G4endl;
#endif

  RotateToLattice(tempvec()=P);
#ifdef G4CMP_DEBUG
  if (verboseLevel>1) G4cout << " in lattice frame " << tempvec() << G4endl;
#endif

  tempvec() = fLattice->MapPToP_Q(ivalley, tempvec());
#ifdef G4CMP_DEBUG
  if (verboseLevel>1) G4cout << " p_q (lattice) " << tempvec() << G4endl;
#endif

  return RotateToSolid(tempvec());
}

G4ThreeVector 
G4LatticePhysical::MapP_QToP(G4int ivalley, const G4ThreeVector& P_Q) const {
#ifdef G4CMP_DEBUG
  if (verboseLevel>1)
    G4cout << "G4LatticePhysical::MapP_QToP " << ivalley << " " << P_Q
	   << G4endl;
#endif

  RotateToLattice(tempvec()=P_Q);
#ifdef G4CMP_DEBUG
  if (verboseLevel>1) G4cout << " in lattice frame " << tempvec() << G4endl;
#endif

  tempvec() = fLattice->MapP_QToP(ivalley, tempvec());
#ifdef G4CMP_DEBUG
  if (verboseLevel>1) G4cout << " p (lattice) " << tempvec() << G4endl;
#endif

  return RotateToSolid(tempvec());
}

G4ThreeVector
G4LatticePhysical::MapV_elToK(G4int ivalley, const G4ThreeVector& v_e) const {
#ifdef G4CMP_DEBUG
  if (verboseLevel>1)
    G4cout << "G4LatticePhysical::MapV_elToK " << ivalley << " " << v_e
     << G4endl;
#endif

  RotateToLattice(tempvec()=v_e);
#ifdef G4CMP_DEBUG
  if (verboseLevel>1) G4cout << " in lattice frame " << tempvec() << G4endl;
#endif

  tempvec() = fLattice->MapV_elToK(ivalley, tempvec());
#ifdef G4CMP_DEBUG
  if (verboseLevel>1) G4cout << " K (lattice) " << tempvec() << G4endl;
#endif

  return RotateToSolid(tempvec());
}

G4ThreeVector 
G4LatticePhysical::MapPtoK(G4int ivalley, const G4ThreeVector& p_e) const {
#ifdef G4CMP_DEBUG
  if (verboseLevel>1)
    G4cout << "G4LatticePhysical::MapPtoK " << ivalley << " " << p_e
     << G4endl;
#endif
  
  RotateToLattice(tempvec()=p_e);
#ifdef G4CMP_DEBUG
  if (verboseLevel>1) G4cout << " in lattice frame " << tempvec() << G4endl;
#endif
  
  tempvec() = fLattice->MapPtoK(ivalley, tempvec());
#ifdef G4CMP_DEBUG
  if (verboseLevel>1) G4cout << " k (lattice) " << tempvec() << G4endl;
#endif

  return RotateToSolid(tempvec());
}

G4ThreeVector 
G4LatticePhysical::MapKtoP(G4int ivalley, const G4ThreeVector& k) const {
#ifdef G4CMP_DEBUG
  if (verboseLevel>1)
    G4cout << "G4LatticePhysical::MapKtoP " << ivalley << " " << k
     << G4endl;
#endif
    
  RotateToLattice(tempvec()=k);
#ifdef G4CMP_DEBUG
  if (verboseLevel>1) G4cout << " in lattice frame " << tempvec() << G4endl;
#endif
    
  tempvec() = fLattice->MapKtoP(ivalley, tempvec());
#ifdef G4CMP_DEBUG
  if (verboseLevel>1) G4cout << " p (lattice) " << tempvec() << G4endl;
#endif

  return RotateToSolid(tempvec());
}

G4double 
G4LatticePhysical::GetElectronEffectiveMass(G4int iv,
					   const G4ThreeVector& p) const {
#ifdef G4CMP_DEBUG
  if (verboseLevel>1)
    G4cout << "G4LatticePhysical::GetElectronEffectiveMass " << iv
	   << " " << p << G4endl;
#endif

  RotateToLattice(tempvec()=p);
#ifdef G4CMP_DEBUG
  if (verboseLevel>1) G4cout << " in lattice frame " << tempvec() << G4endl;
#endif

  return fLattice->GetElectronEffectiveMass(iv, tempvec());
}

G4ThreeVector
G4LatticePhysical::RotateToValley(G4int iv, const G4ThreeVector& v) const {
#ifdef G4CMP_DEBUG
  if (verboseLevel>1)
    G4cout << "G4LatticePhysical::RotateToValley " << iv
	   << " " << v << G4endl;
#endif
  
    RotateToLattice(tempvec()=v);
#ifdef G4CMP_DEBUG
  if (verboseLevel>1) G4cout << " in lattice frame " << tempvec() << G4endl;
#endif
    
    return fLattice->RotateToValley(iv, tempvec());
  }

G4ThreeVector
G4LatticePhysical::RotateFromValley(G4int iv, const G4ThreeVector& v) const {
#ifdef G4CMP_DEBUG
  if (verboseLevel>1)
    G4cout << "G4LatticePhysical::RotateFromValley " << iv
	   << " " << v << G4endl;
#endif

  tempvec() = fLattice->RotateFromValley(iv, v);
  return RotateToSolid(tempvec());
}

G4ThreeVector G4LatticePhysical::
EllipsoidalToSphericalTranformation(G4int iv, const G4ThreeVector& v) const {
#ifdef G4CMP_DEBUG
  if (verboseLevel>1)
    G4cout << "G4LatticePhysical::EllipsoidalToSphericalTranformation " << iv
	   << " " << v << G4endl;
#endif

  RotateToLattice(tempvec()=v);

#ifdef G4CMP_DEBUG
  if (verboseLevel>1) G4cout << " in lattice frame " << tempvec() << G4endl;
#endif

  return fLattice->EllipsoidalToSphericalTranformation(iv, tempvec());
}

// Compute vector in ellipsoidal frame from the spherical frame

G4ThreeVector G4LatticePhysical::
SphericalToEllipsoidalTranformation(G4int iv, const G4ThreeVector& v) const {
#ifdef G4CMP_DEBUG
  if (verboseLevel>1)
    G4cout << "G4LatticePhysical::SphericalToEllipsoidalTranformation " << iv
    << " " << v << G4endl;
#endif

  tempvec() = fLattice->SphericalToEllipsoidalTranformation(iv, v);
  return RotateToSolid(tempvec());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Dump contained logical lattice with volume information

void G4LatticePhysical::Dump(std::ostream& os) const {
  os << "# Physical lattice:"
     << " (hkl) = " << hMiller << " " << kMiller << " " << lMiller
     << " rotation " << fRot/deg << " deg"
     << " @ " << fTemperature/kelvin << " K"
     << "\n# Logical lattice:\n" << *fLattice << std::endl;
}


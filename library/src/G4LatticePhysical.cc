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

#include "G4LatticePhysical.hh"
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
  : verboseLevel(0), fTheta(0), fPhi(0), fLattice(0) {;}

// Set lattice orientation (relative to G4VSolid) with Euler angles

G4LatticePhysical::G4LatticePhysical(const G4LatticeLogical* Lat,
				     G4double theta, G4double phi)
  : verboseLevel(0), fTheta(theta), fPhi(phi), fLattice(Lat) {;}

// Set lattice orientation (relative to G4VSolid) with Miller indices

G4LatticePhysical::G4LatticePhysical(const G4LatticeLogical* Lat,
				     G4int h, G4int k, G4int l)
  : verboseLevel(0), fTheta(0), fPhi(0), fLattice(Lat) {
  SetMillerOrientation(h, k, l);
}

G4LatticePhysical::~G4LatticePhysical() {;}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LatticePhysical::SetLatticeOrientation(G4double t_rot, G4double p_rot) {
  fTheta = t_rot;
  fPhi = p_rot;

  if (verboseLevel) 
    G4cout << "G4LatticePhysical::SetLatticeOrientation " << fTheta << " "
	   << fPhi << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LatticePhysical::SetMillerOrientation(G4int h, G4int k, G4int l) {
  G4ThreeVector norm = h*GetBasis(0) + k*GetBasis(1) + l*GetBasis(2);
  fTheta = norm.theta();
  fPhi = norm.phi();

  if (verboseLevel) 
    G4cout << "G4LatticePhysical::SetMillerOrientation(" << h << k << l 
	   << ") : " << fTheta << " " << fPhi << G4endl;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Rotate input vector between lattice and solid orientations

const G4ThreeVector&
G4LatticePhysical::RotateToLattice(G4ThreeVector& dir) const {
  dir.rotate(CLHEP::HepYHat,fTheta).rotate(CLHEP::HepZHat,fPhi);
  return dir;
}

const G4ThreeVector& 
G4LatticePhysical::RotateToSolid(G4ThreeVector& dir) const {
  dir.rotate(CLHEP::HepZHat,-fPhi).rotate(CLHEP::HepYHat,-fTheta);
  return dir;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

///////////////////////////////
//Loads the group velocity in m/s
/////////////////////////////
G4double G4LatticePhysical::MapKtoV(G4int polarizationState,
				    G4ThreeVector k) const {
  if (verboseLevel>1) G4cout << "G4LatticePhysical::MapKtoV " << k << G4endl;

  RotateToLattice(k);
  return fLattice->MapKtoV(polarizationState, k);
}

///////////////////////////////
//Loads the normalized direction vector along VG
///////////////////////////////
G4ThreeVector G4LatticePhysical::MapKtoVDir(G4int polarizationState,
					    G4ThreeVector k) const {
  if (verboseLevel>1) G4cout << "G4LatticePhysical::MapKtoVDir " << k << G4endl;

  RotateToLattice(k);
  G4ThreeVector VG = fLattice->MapKtoVDir(polarizationState, k);  

  return RotateToSolid(VG);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LatticePhysical::MapPtoEkin(G4int iv, G4ThreeVector p) const {
  if (verboseLevel>1)
    G4cout << "G4LatticePhysical::MapPtoEkin " << iv << " " << p << G4endl;

  RotateToLattice(p);
  return fLattice->MapPtoEkin(iv, p);
}

G4double G4LatticePhysical::MapV_elToEkin(G4int iv, G4ThreeVector v) const {
  if (verboseLevel>1)
    G4cout << "G4LatticePhysical::MapV_elToEkin " << iv << " " << v << G4endl;

  RotateToLattice(v);
  return fLattice->MapV_elToEkin(iv, v);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Convert electron momentum to valley velocity, wavevector, and HV vector

G4ThreeVector 
G4LatticePhysical::MapPtoV_el(G4int ivalley, G4ThreeVector p_e) const {
  if (verboseLevel>1)
    G4cout << "G4LatticePhysical::MapPtoV_el " << ivalley << " " << p_e
	   << G4endl;

  RotateToLattice(p_e);
  p_e = fLattice->MapPtoV_el(ivalley, p_e);	// Overwrite to avoid temporary
  return RotateToSolid(p_e);
}

G4ThreeVector 
G4LatticePhysical::MapV_elToP(G4int ivalley, G4ThreeVector v_e) const {
  if (verboseLevel>1)
    G4cout << "G4LatticePhysical::MapV_elRoP " << ivalley << " " << v_e
	   << G4endl;

  RotateToLattice(v_e);
  v_e = fLattice->MapV_elToP(ivalley, v_e);	// Overwrite to avoid temporary
  return RotateToSolid(v_e);
}

G4ThreeVector
G4LatticePhysical::MapV_elToK_HV(G4int ivalley, G4ThreeVector v_e) const {
  if (verboseLevel>1)
    G4cout << "G4LatticePhysical::MapV_elToK_HV " << ivalley << " " << v_e
     << G4endl;

  RotateToLattice(v_e);
  v_e = fLattice->MapV_elToK_HV(ivalley, v_e);	// Overwrite to avoid temporary
  return RotateToSolid(v_e);
}

G4ThreeVector 
G4LatticePhysical::MapPtoK_valley(G4int ivalley, G4ThreeVector p_e) const {
  if (verboseLevel>1)
    G4cout << "G4LatticePhysical::MapPtoK " << ivalley << " " << p_e
	   << G4endl;

  RotateToLattice(p_e);
  p_e = fLattice->MapPtoK_valley(ivalley, p_e);	// Overwrite to avoid temporary
  return RotateToSolid(p_e);
}

G4ThreeVector 
G4LatticePhysical::MapPtoK_HV(G4int ivalley, G4ThreeVector p_e) const {
  if (verboseLevel>1)
    G4cout << "G4LatticePhysical::MapPtoK_HV " << ivalley << " " << p_e
	   << G4endl;

  RotateToLattice(p_e);
  p_e = fLattice->MapPtoK_HV(ivalley, p_e);	// Overwrite to avoid temporary
  return RotateToSolid(p_e);
}

G4ThreeVector 
G4LatticePhysical::MapK_HVtoK_valley(G4int ivalley, G4ThreeVector k_HV) const {
  if (verboseLevel>1)
    G4cout << "G4LatticePhysical::MapK_HVtoP " << ivalley << " " << k_HV
	   << G4endl;

  RotateToLattice(k_HV);
  k_HV = fLattice->MapK_HVtoK_valley(ivalley, k_HV);
  return RotateToSolid(k_HV);
}

G4ThreeVector
G4LatticePhysical::MapK_HVtoK(G4int ivalley, G4ThreeVector k_HV) const {
  if (verboseLevel>1)
    G4cout << "G4LatticePhysical::MapK_HVtoK " << ivalley << " " << k_HV
     << G4endl;

  RotateToLattice(k_HV);
  k_HV = fLattice->MapK_HVtoK(ivalley, k_HV);	// Overwrite to avoid temporary
  return RotateToSolid(k_HV);
}


G4ThreeVector 
G4LatticePhysical::MapK_HVtoP(G4int ivalley, G4ThreeVector k_HV) const {
  if (verboseLevel>1)
    G4cout << "G4LatticePhysical::MapK_HVtoP " << ivalley << " " << k_HV
	   << G4endl;

  RotateToLattice(k_HV);
  k_HV = fLattice->MapK_HVtoP(ivalley, k_HV);	// Overwrite to avoid temporary
  return RotateToSolid(k_HV);
}

G4ThreeVector 
G4LatticePhysical::MapK_valleyToP(G4int ivalley, G4ThreeVector k) const {
  if (verboseLevel>1)
    G4cout << "G4LatticePhysical::MapK_valleyToP " << ivalley << " " << k
	   << G4endl;

  RotateToLattice(k);
  k = fLattice->MapK_valleyToP(ivalley, k);	// Overwrite to avoid temporary
  return RotateToSolid(k);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Dump contained logical lattice with volume information

void G4LatticePhysical::Dump(std::ostream& os) const {
  os << "# Physical lattice (theta,phi) = "
     << fTheta/deg << " " << fPhi/deg << " deg\n"
     << "# Logical lattice:\n" << *fLattice << std::endl;
}


//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
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

#include "G4LatticePhysical.hh"
#include "G4LatticeLogical.hh"
#include "G4PhysicalConstants.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"


// Unit vectors defined for convenience (avoid memory churn)

namespace {
  G4ThreeVector xhat(1,0,0), yhat(0,1,0), zhat(0,0,1), nullVec(0,0,0);
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
  fTheta = halfpi - std::atan2(l+0.000001,h+0.000001);
  fPhi = halfpi - std::atan2(h+0.000001,k+0.000001);

  if (verboseLevel) 
    G4cout << "G4LatticePhysical::SetMillerOrientation(" << h << k << l 
	   << ") : " << fTheta << " " << fPhi << G4endl;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Rotate input vector between lattice and solid orientations

const G4ThreeVector&
G4LatticePhysical::RotateToLattice(G4ThreeVector& dir) const {
  dir.rotate(yhat,fTheta).rotate(zhat,fPhi);
  return dir;
}

const G4ThreeVector& 
G4LatticePhysical::RotateToSolid(G4ThreeVector& dir) const {
  dir.rotate(zhat,-fPhi).rotate(yhat,-fTheta);
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

// Convert electron momentum to valley velocity

G4ThreeVector 
G4LatticePhysical::MapPtoV_el(G4int ivalley, G4ThreeVector p_e) const {
  if (verboseLevel>1)
    G4cout << "G4LatticePhysical::MapPtoV_el " << ivalley << " " << p_e
	   << G4endl;

  RotateToLattice(p_e);
  const G4RotationMatrix& vToN = GetValley(ivalley);
  G4ThreeVector V = (vToN.inverse()*GetMInvTensor()*vToN) * (p_e/c_light);
  return RotateToSolid(V);
}

// Convert electron momentum to HV wavevector

G4ThreeVector 
G4LatticePhysical::MapPtoK_HV(G4int ivalley, G4ThreeVector p_e) const {
  if (verboseLevel>1)
    G4cout << "G4LatticePhysical::MapPtoK_HV " << ivalley << " " << p_e
	   << G4endl;

  RotateToLattice(p_e);
  p_e.transform(GetValley(ivalley));		// Rotate into valley frame
  G4ThreeVector k_HV = GetSqrtInvTensor() * (p_e/hbarc);    // Henning-Vogt
  return RotateToSolid(k_HV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Dump contained logical lattice with volume information

void G4LatticePhysical::Dump(std::ostream& os) const {
  os << "# Physical lattice (theta,phi) = "
     << fTheta/deg << " " << fPhi/deg << " deg\n"
     << "# Logical lattice:\n" << *fLattice << std::endl;
}


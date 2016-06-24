//  G4CMPPhononKinematics.cpp
//  Created by Daniel Palken in 2014 for G4CMP
//
//  20160624  Allow non-unit vector to be passed into computeKinematics()

#include "G4CMPPhononKinematics.hh"
#include "G4LatticeLogical.hh"
#include "G4PhononPolarization.hh"
#include "G4ThreeVector.hh"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
using namespace std;

// ++++++++++++++++++++++ G4CMPPhononKinematics METHODS +++++++++++++++++++++++++++

G4CMPPhononKinematics::G4CMPPhononKinematics(G4LatticeLogical *lat)
  : lattice(lat), christoffel(G4ThreeVector::SIZE, G4ThreeVector::SIZE, 0.) {;}

G4CMPPhononKinematics::~G4CMPPhononKinematics() {;}

// Build D_il, the Christoffel matrix that defines the eigensystem
void G4CMPPhononKinematics::fillChristoffelMatrix(const G4ThreeVector& nn)
{
  christoffel.clear();
  for (int i = 0; i < G4ThreeVector::SIZE; i++) {
    for (int l = 0; l < G4ThreeVector::SIZE; l++) {
      for (int j = 0; j < G4ThreeVector::SIZE; j++) {
	for (int m = 0; m < G4ThreeVector::SIZE; m++) {
	  christoffel[i][l] += (lattice->GetCijkl(i,j,l,m) * nn[j] * nn[m]);
	}
      }
      christoffel[i][l] /= lattice->GetDensity();
    }
  }
}

// Compute kinematics for specified wavevector (direction)
void G4CMPPhononKinematics::computeKinematics(const G4ThreeVector& n_dir) {
  if (!n_dir.unit().isNear(last_ndir)) return;		// Already computed

  /* get the Christoffel Matrix D_il, which is symmetric (it
     equals its transpose).  This also means its eigenvalues will
     all be real (NR, pg. 564) */
  fillChristoffelMatrix(n_dir.unit());
  
  /* set up and solve eigensystem of D_il:
     Use NR's method for real, symmetric matricies.
     Eigenvalues are the phase velocities squared (v_phase = omega/k).
     Eigenvectors are the corresponding polaizrations e_l.
     Eigenvalues stored in eigenSys.d[0..n-1] in descening order.
     Corresponding eigenvectors are the columns of eigenSys.z[0..n-1][0..n-1] */
  eigenSys.setup(christoffel);
  
  /* Extract eigen vectors and values for each mode */
  for (int mode = 0; mode < G4PhononPolarization::NUM_MODES; mode++) {
    // calculate desired quantities that will populate lookup table:
    vphase[mode] = sqrt(eigenSys.d[mode]);
    slowness[mode] = n_dir.unit()/vphase[mode];
    polarization[mode].set(eigenSys.z[G4ThreeVector::X][mode],
			   eigenSys.z[G4ThreeVector::Y][mode],
			   eigenSys.z[G4ThreeVector::Z][mode]);
    
    computeGroupVelocity(mode, eigenSys.z, slowness[mode]);
  }
  
  /* Store wavevector direction to avoid recalculations */
  last_ndir = n_dir.unit();
}

// Fill group velocity cache for specified mode from lattice parameters
// NOTE:  Must only be called from computeKinematics() above!
void G4CMPPhononKinematics::computeGroupVelocity(int mode, const matrix<double>& e_mat,
					     const G4ThreeVector& slow) {
  vgroup[mode].set(0.,0.,0.);
  for (int dim=0; dim<G4ThreeVector::SIZE; dim++) {
    for (int i=0; i<G4ThreeVector::SIZE; i++) {
      for (int j=0; j<G4ThreeVector::SIZE; j++) {
	for (int l=0; l<G4ThreeVector::SIZE; l++) {
	  vgroup[mode][dim] += (e_mat[i][mode] * lattice->GetCijkl(i,j,l,dim)
				* slow[j] * e_mat[l][mode]);
	}
      }
    }
  }
  
  vgroup[mode] /= lattice->GetDensity();
}

const G4ThreeVector& 
G4CMPPhononKinematics::getGroupVelocity(int mode, const G4ThreeVector& n_dir) {
  computeKinematics(n_dir);
  return vgroup[mode];
}

const G4ThreeVector& 
G4CMPPhononKinematics::getPolarization(int mode, const G4ThreeVector& n_dir) {
  computeKinematics(n_dir);
  return polarization[mode];
}

const G4ThreeVector& 
G4CMPPhononKinematics::getSlowness(int mode, const G4ThreeVector& n_dir) {
  computeKinematics(n_dir);
  return slowness[mode];
}

double 
G4CMPPhononKinematics::getPhaseSpeed(int mode, const G4ThreeVector& n_dir) {
  computeKinematics(n_dir);
  return vphase[mode];
}

const G4String& G4CMPPhononKinematics::getLatticeName() const {
  return lattice->GetName();
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

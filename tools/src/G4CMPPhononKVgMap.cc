//  G4CMPPhononKVgMap.cpp
//  Created by Daniel Palken in 2014 for G4CMP

#include "G4CMPPhononKVgMap.hh"
#include "G4LatticeLogical.hh"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
using namespace std;

// '''''''''''''''''''''''''''''' PRIVATE CONSTANTS ''''''''''''''''''''''''''''

// other
#define REDUCED_TENSOR_SIZE     6           // 2D tensor, 36 elements
#define FULL_TENSOR_SIZE        3           // 4D tensor, 81 elements

// cartesian indexing
#define X_DIR                   0
#define Y_DIR                   1
#define Z_DIR                   2
// in total:
#define SPATIAL_DIMENSIONS      3           // x, y, z

// ++++++++++++++++++++++ G4CMPPhononKVgMap STRUCT METHODS +++++++++++++++++++++

G4CMPPhononKVgMap::G4CMPPhononKVgMap(G4LatticeLogical *lat)
  : lattice(lat), christoffel(SPATIAL_DIMENSIONS, SPATIAL_DIMENSIONS, 0.) {;}

G4CMPPhononKVgMap::~G4CMPPhononKVgMap() {;}


// ****************************** BUILD METHODS ********************************

// Build D_il, the Christoffel matrix that defines the eigensystem
void G4CMPPhononKVgMap::fillChristoffelMatrix(const G4ThreeVector& nn)
{
  christoffel.clear();
  for (int i = 0; i < SPATIAL_DIMENSIONS; i++) {
    for (int l = 0; l < SPATIAL_DIMENSIONS; l++) {
      for (int j = 0; j < SPATIAL_DIMENSIONS; j++) {
	for (int m = 0; m < SPATIAL_DIMENSIONS; m++) {
	  christoffel[i][l] += (lattice->GetCijkl(i,j,l,m) * nn[j] * nn[m]);
	}
      }
      christoffel[i][l] /= lattice->GetDensity();
    }
  }
}

// Compute kinematics for specified wavevector (direction)
void G4CMPPhononKVgMap::computeKinematics(const G4ThreeVector& n_dir) {
  /* get the Christoffel Matrix D_il, which is symmetric (it
     equals its transpose).  This also means its eigenvalues will
     all be real (NR, pg. 564) */
  fillChristoffelMatrix(n_dir);
  
  /* set up and solve eigensystem of D_il:
     Use NR's method for real, symmetric matricies.
     Eigenvalues are the phase velocities squared (v_phase = omega/k).
     Eigenvectors are the corresponding polaizrations e_l.
     Eigenvalues stored in eigenSys.d[0..n-1] in descening order.
     Corresponding eigenvectors are the columns of eigenSys.z[0..n-1][0..n-1] */
  eigenSys.setup(christoffel);
  
  /* Extract eigen vectors and values for each mode */
  for (int mode = 0; mode < NUM_MODES; mode++) {
    // calculate desired quantities that will populate lookup table:
    vphase[mode] = sqrt(eigenSys.d[mode]);
    slowness[mode] = n_dir/vphase[mode];
    polarization[mode].set(eigenSys.z[X_DIR][mode], eigenSys.z[Y_DIR][mode],
			   eigenSys.z[Z_DIR][mode]);
    
    computeGroupVelocity(mode, eigenSys.z, slowness[mode]);
  }
  
  /* Store wavevector direction to avoid recalculations */
  last_ndir = n_dir;
}

// Fill group velocity cache for specified mode from lattice parameters
void G4CMPPhononKVgMap::computeGroupVelocity(int mode, const matrix<double>& e_mat,
					     const G4ThreeVector& slow) {
  vgroup[mode].set(0.,0.,0.);
  for (int dim=0; dim<SPATIAL_DIMENSIONS; dim++) {
    for (int i=0; i<FULL_TENSOR_SIZE; i++) {
      for (int j=0; j<FULL_TENSOR_SIZE; j++) {
	for (int l=0; l<FULL_TENSOR_SIZE; l++) {
	  vgroup[mode][dim] += (e_mat[i][mode] * lattice->GetCijkl(i,j,l,dim)
				* slow[j] * e_mat[l][mode]);
	}
      }
    }
  }
  
  vgroup[mode] /= lattice->GetDensity();
}

const G4ThreeVector& G4CMPPhononKVgMap::getGroupVelocity(int mode, const G4ThreeVector& n_dir) {
  if (!n_dir.isNear(last_ndir)) computeKinematics(n_dir);
  return vgroup[mode];
}

const G4ThreeVector& G4CMPPhononKVgMap::getPolarization(int mode, const G4ThreeVector& n_dir) {
  if (!n_dir.isNear(last_ndir)) computeKinematics(n_dir);
  return polarization[mode];
}

const G4ThreeVector& G4CMPPhononKVgMap::getSlowness(int mode, const G4ThreeVector& n_dir) {
  if (!n_dir.isNear(last_ndir)) computeKinematics(n_dir);
  return slowness[mode];
}

double G4CMPPhononKVgMap::getPhaseSpeed(int mode, const G4ThreeVector& n_dir) {
  if (!n_dir.isNear(last_ndir)) computeKinematics(n_dir);
  return vphase[mode];
}

const G4String& G4CMPPhononKVgMap::getLatticeName() const {
  return lattice->GetName();
}

// *****************************************************************************

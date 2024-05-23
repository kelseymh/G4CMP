//  G4CMPPhononKinematics.cpp
//  Created by Daniel Palken in 2014 for G4CMP
//
//  20160624  Allow non-unit vector to be passed into computeKinematics()
//  20170525  Drop unnecessary empty destructor ("rule of five" semantics)

#include "G4CMPPhononKinematics.hh"
#include "G4LatticeLogical.hh"
#include "G4PhononPolarization.hh"
#include "G4ThreeVector.hh"

// ++++++++++++++++++++++ G4CMPPhononKinematics METHODS +++++++++++++++++++++++++++

G4CMPPhononKinematics::G4CMPPhononKinematics(G4LatticeLogical *lat)
  : lattice(lat), christoffel(G4ThreeVector::SIZE, G4ThreeVector::SIZE, 0.) {;}

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
  if (n_dir.unit().isNear(last_ndir)) return;		// Already computed

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
  
  /* Extract eigen vectors and values for each mode.
   * We must sort them to match the sorting in G4PhononPolarization.
   * This assumes that fast transverse is more energetic than slow transverse,
   * and I'm not positive that's always true.
   */

  //Whichever eigenvector is most parallel with k is the longitudinal mode
  G4double mostParallelMeasure = 0;
  size_t longIdx = 0;
  for (size_t i = 0; i < 3; ++i) {
    const G4double howParallel = G4ThreeVector(eigenSys.z[0][i],
                                               eigenSys.z[1][i],
                                               eigenSys.z[2][i])
                                              .howOrthogonal(n_dir);
    if (howParallel > mostParallelMeasure) {
      mostParallelMeasure = howParallel;
      longIdx = i;
    }
  }

  // For the tranverse, we assume the smallest eigenvalue is slow transverse
  // Lucky for us (kinda), the eigenvalues come out sorted.
  size_t slowTransIdx = 0;
  size_t fastTransIdx = 0;
  for (size_t i = 0; i < 3; ++i) {
    if (i != longIdx) {
      fastTransIdx = i;
      slowTransIdx = (i + 1) % 3 == longIdx ? (i + 2) % 3 : (i + 1) % 3;
      break;
    }
  }

  for (int mode = 0; mode < G4PhononPolarization::NUM_MODES; mode++) {
    // calculate desired quantities that will populate lookup table.
    // Must map the G4PhononPolarization indices to the eigen indices from
    // the solver.
    size_t idx = (mode == G4PhononPolarization::Long ? longIdx :
		  mode == G4PhononPolarization::TransFast ? fastTransIdx :
		  slowTransIdx);
    vphase[mode] = sqrt(eigenSys.d[idx]);
    slowness[mode] = n_dir.unit()/vphase[mode];
    polarization[mode].set(eigenSys.z[G4ThreeVector::X][idx],
                           eigenSys.z[G4ThreeVector::Y][idx],
                           eigenSys.z[G4ThreeVector::Z][idx]);
    
    computeGroupVelocity(mode, idx, eigenSys.z, slowness[mode]);
  }
  
  /* Store wavevector direction to avoid recalculations */
  last_ndir = n_dir.unit();
}

// Fill group velocity cache for specified mode from lattice parameters
// NOTE:  Must only be called from computeKinematics() above!
void G4CMPPhononKinematics::computeGroupVelocity(int mode,
                                                 size_t idx,
                                                 const matrix<double>& e_mat,
                                                 const G4ThreeVector& slow) {
  vgroup[mode].set(0.,0.,0.);
  for (int dim=0; dim<G4ThreeVector::SIZE; dim++) {
    for (int i=0; i<G4ThreeVector::SIZE; i++) {
      for (int j=0; j<G4ThreeVector::SIZE; j++) {
	for (int l=0; l<G4ThreeVector::SIZE; l++) {
    vgroup[mode][dim] += (e_mat[i][idx] * lattice->GetCijkl(i,j,l,dim)
        * slow[j] * e_mat[l][idx]);
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

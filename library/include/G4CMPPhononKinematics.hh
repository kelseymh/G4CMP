#ifndef G4CMPPhononKinematics_hh
#define G4CMPPhononKinematics_hh

//  G4CMPPhononKinematics.hh
//  Created by Daniel Palken in 2014 for G4CMP
//
//  20170525  Drop unnecessary empty destructor ("rule of five" semantics)

#include "G4CMPEigenSolver.hh" // Numerical Recipes III code
#include "G4CMPMatrix.hh"
#include "G4PhononPolarization.hh"
#include "G4ThreeVector.hh"
#include <string>
#include <vector>
using std::string;
using std::vector;
using G4CMP::matrix;

class G4LatticeLogical;


class G4CMPPhononKinematics {
public:
  G4CMPPhononKinematics(G4LatticeLogical *lat);

  // Direct calculations
  void computeKinematics(const G4ThreeVector& n_dir);
  void fillChristoffelMatrix(const G4ThreeVector& n_dir);
  void computeGroupVelocity(int mode, size_t idx, const matrix<double>& epol,
			    const G4ThreeVector& slow);
  const G4ThreeVector& getGroupVelocity(int mode, const G4ThreeVector& n_dir);
  const G4ThreeVector& getPolarization(int mode, const G4ThreeVector& n_dir);
  const G4ThreeVector& getSlowness(int mode, const G4ThreeVector& n_dir);
  double getPhaseSpeed(int mode, const G4ThreeVector& n_dir);

public:
  const G4String& getLatticeName() const;	// For use with lookup table

private:
  G4LatticeLogical* lattice;

  // Data buffers to compute kinematics for all modes in specified direction
  G4ThreeVector last_ndir;		// Buffer to handle caching results
  G4CMPEigenSolver eigenSys;
  matrix<double> christoffel;
  double vphase[G4PhononPolarization::NUM_MODES];
  G4ThreeVector slowness[G4PhononPolarization::NUM_MODES];
  G4ThreeVector vgroup[G4PhononPolarization::NUM_MODES];
  G4ThreeVector polarization[G4PhononPolarization::NUM_MODES];
};

#endif /* G4CMPPhononKinematics_hh */

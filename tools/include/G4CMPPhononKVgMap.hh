//  G4CMPPhononKVgMap.hh
//  Created by Daniel Palken in 2014 for G4CMP

#ifndef G4CMPPhononKVgMap_hh
#define G4CMPPhononKVgMap_hh

#include "G4CMPEigenSolver.hh" // Numerical Recipes III code
#include "G4CMPInterpolator.hh"
#include "G4LatticeLogical.hh"
#include "G4ThreeVector.hh"
#include "matrix.hh"
#include <string>
#include <vector>
using std::string;
using std::vector;
using G4CMP::matrix;

// '''''''''''''''''''''''''''''' PUBLIC CONSTANTS ''''''''''''''''''''''''''''''''

// magic values for use with interpolation
#define OUT_OF_BOUNDS           9.0e299     // signals looking outside of the lookup table
#define ERRONEOUS_INPUT         1.0e99      // signals that input is not correct
// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

// ++++++++++++++++++++++++++++++ G4CMPPhononKVgMap STRUCT ++++++++++++++++++++++++
struct G4CMPPhononKVgMap {
public:
  // Wave mode indexing (longitudinal, slow and fast transverse)
  enum PhononModes { L, ST, FT, NUM_MODES };
  string getModeName(int MODE);
    
public:
  G4CMPPhononKVgMap(G4LatticeLogical *lat);
  ~G4CMPPhononKVgMap();

  // Direct calculations
  void computeKinematics(const G4ThreeVector& n_dir);
  void fillChristoffelMatrix(const G4ThreeVector& n_dir);
  void computeGroupVelocity(int mode, const matrix<double>& epol,
			    const G4ThreeVector& slow);
  const G4ThreeVector& getGroupVelocity(int mode, const G4ThreeVector& n_dir);
  const G4ThreeVector& getPolarization(int mode, const G4ThreeVector& n_dir);
  const G4ThreeVector& getSlowness(int mode, const G4ThreeVector& n_dir);
  double getPhaseSpeed(int mode, const G4ThreeVector& n_dir);

public:
  // Symbolic identifiers for various arrays
  enum DataTypes { N_X, N_Y, N_Z, THETA, PHI,	// Wavevector
		   S_X, S_Y, S_Z, S_MAG, S_PAR,	// Slowness (1/vphase)
		   V_P,				// Vphase = omega/k
		   V_G, V_GX, V_GY, V_GZ,	// Vgroup
		   E_X, E_Y, E_Z,		// Polarization
		   NUM_DATA_TYPES };
  string getDataTypeName(int TYPE);

  // interpolation methods
  double interpolateEven(double nx, double ny, int MODE, int TYPE_OUT, bool SILENT=true);
  double interpolateEven(G4CMPBiLinearInterp& grid, double nx, double ny);
  double interpGeneral(int mode, const G4ThreeVector& k, int typeDesired);

  G4ThreeVector interpGroupVelocity_N(int mode, const G4ThreeVector& k);

  double interpPerpSlowness(int mode, const G4ThreeVector& k)
  { return interpGeneral(mode, k, S_Z); }

  double interpGroupVelocity(int mode, const G4ThreeVector& k)
  { return interpGeneral(mode, k, V_G); }
  
  // Dump lookup table for external use
  void writeLookupTable();

private:
  // Populate full table for interpolation
  void setUpDataVectors();
  void generateLookupTable();
  void generateMultiEvenTable();
  G4CMPBiLinearInterp generateEvenTable(PhononModes MODE, DataTypes TYPE_OUT);
  void clearQuantityMap();

private:
  G4LatticeLogical* lattice;
  vector<vector<G4CMPBiLinearInterp> > quantityMap;
  vector<vector<vector<double> > > lookupData;

private:
  // Data buffers to compute kinematics for all modes in specified direction
  G4ThreeVector last_ndir;		// Buffer to handle caching results
  G4CMPEigenSolver eigenSys;
  matrix<double> christoffel;
  double  vphase[NUM_MODES];
  G4ThreeVector slowness[NUM_MODES];
  G4ThreeVector vgroup[NUM_MODES];
  G4ThreeVector polarization[NUM_MODES];
};
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// ############################# UTILITY FUNCTIONS ################################
bool doubApproxEquals(double D1, double D2, double Tol = 1.0e-10,
                      bool USE_ABS_TOL = true);
bool doubGreaterThanApprox(double D1, double D2, double TOL = 1.0e-10,
                           bool USE_ABS_TOL = true);
bool doubLessThanApprox(double D1, double D2, double TOL = 1.0e-10,
                        bool USE_ABS_TOL = true);
// #################################################################################

#endif /* G4CMPPhononKVgMap_hh */

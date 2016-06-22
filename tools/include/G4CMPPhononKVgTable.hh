//  G4CMPPhononKVgTable.hh
//  Created by Daniel Palken in 2014 for G4CMP

#ifndef G4CMPPhononKVgTable_hh
#define G4CMPPhononKVgTable_hh

#include "G4CMPInterpolator.hh"
#include "G4ThreeVector.hh"
#include <string>
#include <vector>
using std::string;
using std::vector;

class G4CMPPhononKVgMap;


// '''''''''''''''''''''''''''''' PUBLIC CONSTANTS ''''''''''''''''''''''''''''''''

// magic values for use with interpolation
#define OUT_OF_BOUNDS           9.0e299     // signals looking outside of the lookup table
#define ERRONEOUS_INPUT         1.0e99      // signals that input is not correct
// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

// ++++++++++++++++++++++++++++++ G4CMPPhononKVgTable ++++++++++++++++++++++++
class G4CMPPhononKVgTable {
public:
  G4CMPPhononKVgTable(G4CMPPhononKVgMap* map);
  ~G4CMPPhononKVgTable();

public:
  // Symbolic identifiers for various arrays, to use with lookup table
  enum DataTypes { N_X, N_Y, N_Z, THETA, PHI,	// Wavevector
		   S_X, S_Y, S_Z, S_MAG, S_PAR,	// Slowness (1/vphase)
		   V_P,				// Vphase = omega/k
		   V_G, V_GX, V_GY, V_GZ,	// Vgroup
		   E_X, E_Y, E_Z,		// Polarization
		   NUM_DATA_TYPES };
  string getDataTypeName(int TYPE);
  string getModeName(int MODE);

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
  void write();

private:
  // Populate full table for interpolation
  void setUpDataVectors();
  void generateLookupTable();
  void generateMultiEvenTable();
  G4CMPBiLinearInterp generateEvenTable(int MODE, DataTypes TYPE_OUT);
  void clearQuantityMap();

private:
  G4CMPPhononKVgMap* mapper;
  vector<vector<G4CMPBiLinearInterp> > quantityMap;
  vector<vector<vector<double> > > lookupData;
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

#endif /* G4CMPPhononKVgTable_hh */

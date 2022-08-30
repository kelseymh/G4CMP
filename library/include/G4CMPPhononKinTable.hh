//  G4CMPPhononKinTable.hh
//  Created by Daniel Palken in 2014 for G4CMP
//
//  20160628  Tabulating on nx and ny is just wrong; use theta, phi
//  20170525  Drop unnecessary empty destructor ("rule of five" semantics)

#ifndef G4CMPPhononKinTable_hh
#define G4CMPPhononKinTable_hh

#include "G4CMPInterpolator.hh"
#include "G4PhysicalConstants.hh"
#include "G4ThreeVector.hh"
#include <string>
#include <vector>
using std::string;
using std::vector;

class G4CMPPhononKinematics;


// ++++++++++++++++++++++++++++++ G4CMPPhononKinTable +++++++++++++++++++++++++
class G4CMPPhononKinTable {
public:
  G4CMPPhononKinTable(G4CMPPhononKinematics* map, G4double thmin=0.,
		      G4double thmax=pi, G4int nth=250, G4double phmin=0.,
		      G4double phmax=twopi, G4int nph=250);

  void initialize();		// Trigger filling of lookup tables

public:
  // Symbolic identifiers for various arrays, to use with lookup table
  enum DataTypes { N_X, N_Y, N_Z, THETA, PHI,	// Wavevector
		   S_X, S_Y, S_Z, S_MAG, S_PAR,	// Slowness (1/vphase)
		   V_P,				// Vphase = omega/k
		   V_G, V_GX, V_GY, V_GZ,	// Vgroup
		   E_X, E_Y, E_Z,		// Polarization
		   NUM_DATA_TYPES };
  string getDataTypeName(int TYPE);
  G4double getDataUnit(int TYPE);	   // Geant4 units for writing output

  // Special values for interpolation errors
  static const G4double OUT_OF_BOUNDS;     // Looking outside of lookup table
  static const G4double ERRONEOUS_INPUT;   // Input is not correct

  bool goodBin(G4double theta, G4double phi) {
    return (thetaMin <= theta && theta <= thetaMax &&
	    phiMin <= phi && phi <= phiMax);
  }

  // interpolation methods
  double interpGeneral(int mode, const G4ThreeVector& k, int typeDesired);

  G4ThreeVector interpGroupVelocity_N(int mode, const G4ThreeVector& k);

  double interpPerpSlowness(int mode, const G4ThreeVector& k)
  { return interpGeneral(mode, k, S_Z); }

  double interpGroupVelocity(int mode, const G4ThreeVector& k)
  { return interpGeneral(mode, k, V_G); }
  
  // Dump lookup table for external use
  void write();

protected:
  // Internal drivers for lookup tables
  double interpolateEven(double theta, double phi, int MODE, int TYPE_OUT,
			 bool SILENT=true);
  double interpolateEven(G4CMPGridInterp& grid, double theta, double phi);

private:
  G4double thetaMin, thetaMax, thetaStep;   // Range and steps for wavevector
  G4int thetaCount;
  G4double phiMin, phiMax, phiStep;
  G4int phiCount;

  // Populate full table for interpolation
  void setUpDataVectors();
  void generateLookupTable();
  void generateMultiEvenTable();
  G4CMPGridInterp generateEvenTable(int MODE, DataTypes TYPE_OUT);
  void clearQuantityMap();

private:
  G4CMPPhononKinematics* mapper;	// Not owned; client responsibility
  G4bool lookupReady;			// Flag once tables are filled
  vector<vector<G4CMPGridInterp> > quantityMap;
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

#endif /* G4CMPPhononKinTable_hh */

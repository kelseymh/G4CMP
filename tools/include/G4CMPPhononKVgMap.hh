//  G4CMPG4CMPPhononKVgMap.h
//  Created by Daniel Palken in 2014 for G4CMP

#ifndef __G4CMPG4CMPPhononKVgMap__
#define __G4CMPG4CMPPhononKVgMap__

#include "G4CMPEigenSolver.hh" // Numerical Recipes III code
#include "G4ThreeVector.hh"
#include <string>
using namespace std;

// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
// '''''''''''''''''''''''''''''' PUBLIC CONSTANTS ''''''''''''''''''''''''''''''''
// data type indexing
/* note that this is based on the index (starting at 0) within the large data structure,
 and is off by one from the column number (still starting at 0) in the lookup table
 itself, due to the mode column which is present to make the lookup table a good deal
 more read- and debugg-able */
#define N_X                     0           // k_x directional component
#define N_Y                     1           // k_y directional component
#define N_Z                     2           // k_z directional component
#define THETA                   3           // angle out from x-axis in xy-plane
#define PHI                     4           // angle down from z-axis
#define S_X                     5           // slowness x component
#define S_Y                     6           // slowness y component
#define S_Z                     7           // slowness z component
#define S_MAG                   8           // slowness magnitude
#define S_PAR                   9           // slowness parallel component
#define V_P                     10          // phase velocity = omega/k
#define V_G                     11          // group velocity = grad_k(omega)
#define V_GX                    12          // group velocity x component
#define V_GY                    13          // group velocity y component
#define V_GZ                    14          // group velocity z component
#define E_X                     15          // polarization x compoent
#define E_Y                     16          // polarization y compoent
#define E_Z                     17          // polarization z compoent
// in total:
#define NUM_DATA_TYPES          18          // see above

// wave mode indexing for eigenvalues and vectors
#define L                       0           // longitudinal
#define ST                      1           // slow transverse
#define FT                      2           // fast transverse
// in total:
#define NUM_MODES               3           // L, ST, FT

// interpolation
#define OUT_OF_BOUNDS           9.0e299     // signals looking outside of the lookup table
#define ERRONEOUS_INPUT         1.0e99      // signals that input is not correct
// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''


// ||||||||||||||||||||||||| Material STRUCT ||||||||||||||||||||||||||||||||||||||
struct Material {
private:
  typedef double**** CMatrix;		// For elasticity matrix below

  string NAME;				// the name of the material
  double DENSITY;			// the density of the material in kg/m^3
  double C_11, C_12, C_44;		// assumes a high level of symmetry
  MatDoub C_reduced;			// Second-order elasticity tensor
  CMatrix C_full;			// Complete fourth-order elasticity tensor
  MatInt rn, rn2;			// Index mapping for tensor generation

public:
  Material(const string& MATNAME, double DEN, double C11, double C12, double C44);
  ~Material();

  // getters:
  const string& getName() {return NAME;}
  double getDensity() {return DENSITY;}
  CMatrix getC_full() {return C_full;}
  double getC_ijlm(int i, int j, int l, int m) {return C_full[i][j][l][m];}
  MatDoub getChristoffelMatrix(const G4ThreeVector& n_dir);
  
private:
  void fillReducedTensor();
  void fillFullTensor();
  void generate_rn();       // Index mapping between full and reduced tensors
  void generate_rn2();
};
// ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++ G4CMPPhononKVgMap STRUCT ++++++++++++++++++++++++++++++++++
struct G4CMPPhononKVgMap {
private:
  Material* material;
  vector<vector<G4CMPBiLinearInterp> > quantityMap;
  vector<vector<vector<double> > > lookupData;
    
public:
  G4CMPPhononKVgMap(Material *mat) : material(mat) {
    generateLookupTable();
    generateMultiEvenTable();
  }

  ~G4CMPPhononKVgMap() {
    clearQuantityMap();
  }

  // Direct calculations
  void computeKinematics(const G4ThreeVector& n_dir);
  void computeGroupVelocity(int mode, const MatDoub& epol,
			    const G4ThreeVector& slow);
  const G4ThreeVector& getGroupVelocity(int mode, const G4ThreeVector& n_dir);
  const G4ThreeVector& getPolarization(int mode, const G4ThreeVector& n_dir);
  const G4ThreeVector& getSlowness(int mode, const G4ThreeVector& n_dir);
  double getPhaseSpeed(int mode, const G4ThreeVector& n_dir);

  // interpolation methods
  double interpolateEven(double nx, double ny, int MODE, int TYPE_OUT, bool SILENT=true);
  double interpolateEven(G4CMPBiLinearInterp& grid, double nx, double ny);
  double interpGeneral(int mode, const G4ThreeVector& k, int typeDesired);

  G4ThreeVector interpGroupVelocity_N(int mode, const G4ThreeVector& k);

  double interpPerpSlowness(int mode, const G4ThreeVector& k)
  { return interpGeneral(mode, k, S_Z); }

  double interpGroupVelocity(int mode, const G4ThreeVector& k)
  { return interpGeneral(mode, k, V_G); }
  
  // getter
  const string& getMatName() { return material->getName(); }

  // Dump lookup table for external use
  void writeLookupTable();

private:
  // Populate full table for interpolation
  void setUpDataVectors();
  void generateLookupTable();
  void generateMultiEvenTable();
  G4CMPBiLinearInterp generateEvenTable(int MODE, int TYPE_OUT);
  void clearQuantityMap();

private:
  // Data buffers to compute kinematics for all modes in specified direction
  G4ThreeVector last_ndir;		// Buffer to handle caching results
  G4CMPEigenSolver eigenSys;
  MatDoub christoffel;
  double  vphase[NUM_MODES];
  G4ThreeVector slowness[NUM_MODES];
  G4ThreeVector vgroup[NUM_MODES];
  G4ThreeVector polarization[NUM_MODES];
};
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// ############################# EXTRANEOUS HEADERS ###############################
string getModeName(int MODE);
string getDataTypeName(int DATA_TYPE);

bool doubApproxEquals(double D1, double D2, double Tol = 1.0e-10,
                      bool USE_ABS_TOL = true);
bool doubGreaterThanApprox(double D1, double D2, double TOL = 1.0e-10,
                           bool USE_ABS_TOL = true);
bool doubLessThanApprox(double D1, double D2, double TOL = 1.0e-10,
                        bool USE_ABS_TOL = true);
// #############################################################################

#endif /* defined(__KMappingCondensed__G4CMPPhononKVgMap__) */

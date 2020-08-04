//  G4CMPPhononKinTable.cpp
//  Created by Daniel Palken in 2014 for G4CMP
//
//  20160627  M. Kelsey -- Defer table building to first query
//  20160628  Tabulating on nx and ny is just wrong; use theta, phi
//  20170525  Drop unnecessary empty destructor ("rule of five" semantics)
//  20170527  Abort job if output file fails

#include "G4CMPPhononKinTable.hh"
#include "G4CMPMatrix.hh"
#include "G4CMPPhononKinematics.hh"
#include "G4PhononPolarization.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace std;
using G4CMP::matrix;

// ++++++++++++++++++++++ G4CMPPhononKinTable METHODS +++++++++++++++++++++++++

// magic values for use with interpolation
const G4double G4CMPPhononKinTable::OUT_OF_BOUNDS = 9.0e299;
const G4double G4CMPPhononKinTable::ERRONEOUS_INPUT = 1.0e99;

G4CMPPhononKinTable::
G4CMPPhononKinTable(G4CMPPhononKinematics* map,
		    G4double thmin, G4double thmax, G4int nth,
		    G4double phmin, G4double phmax, G4int nph)
  : thetaMin(thmin), thetaMax(thmax),
    thetaStep((nth>0)?(thmax-thmin)/nth:1.), thetaCount(nth),
    phiMin(phmin), phiMax(phmax),
    phiStep((nph>0)?(phmax-phmin)/nph:1.), phiCount(nph),
    mapper(map), lookupReady(false) {;}

void G4CMPPhononKinTable::initialize() {
  if (lookupReady) return;		// Tables already generated

  generateLookupTable();
  generateMultiEvenTable();
  lookupReady = true;
}

void G4CMPPhononKinTable::clearQuantityMap() {
  // common technique to free up the memory of a vector:
  quantityMap.clear(); // this alone actually does not free up the memory
  vector<vector<G4CMPGridInterp> > newVec;
  quantityMap.swap(newVec);

  lookupReady = false;
}

// returns aphi quantity desired from the interpolation table
double G4CMPPhononKinTable::interpGeneral(int mode, const G4ThreeVector& k,
					  int typeDesired) {
  if (!lookupReady) initialize();	// Fill tables on first query

  // Angles must be in range [0,pi) and [0,twopi)
  double theta = k.theta(); theta+=(theta<0.)?pi:0.;
  double phi = k.phi();     phi+=(phi<0.)?twopi:0.;

  // note: this method at present does not require nz
  return interpolateEven(theta, phi, mode, typeDesired);
}

// returns the unit vector pointing in the direction of Vg
G4ThreeVector 
G4CMPPhononKinTable::interpGroupVelocity_N(int mode, const G4ThreeVector& k) {
  G4ThreeVector Vg(interpGeneral(mode, k, V_GX),
		   interpGeneral(mode, k, V_GY),
		   interpGeneral(mode, k, V_GZ));
  return Vg.unit();
}

// ****************************** BUILD METHODS ********************************
/* sets up the vector of vectors of vectors used to store the data
   from the lookup table */
void G4CMPPhononKinTable::setUpDataVectors() {
  // Preload outer vectors with correct structure, to avoid push-backs
  lookupData.resize(G4PhononPolarization::NUM_MODES,
		    vector<vector<double> >(NUM_DATA_TYPES));
}

// makes the lookup table for whatever material is specified
void G4CMPPhononKinTable::generateLookupTable() {
  setUpDataVectors();

  // Kinematic data buffers fetched from Mapper (avoids memory churn)
  G4double vphase;
  G4ThreeVector slowness;
  G4ThreeVector vgroup;
  G4ThreeVector polarization;

#ifdef G4CMP_DEBUG
  cout << "G4CMPPhononKinTable: "
       << thetaCount << " X bins [" << thetaMin << ".." << thetaMax << "], "
       << phiCount << " X bins [" << phiMin << ".." << phiMax << "]"
       << G4endl;
#endif

  // ensures even spacing for x and y on a unit circle in the xy-plane
  G4ThreeVector n_dir;
  for (int ith=0; ith<=thetaCount; ith++) {
    double theta = thetaMin + ith*thetaStep;

    for (int iphi=0; iphi<=phiCount; iphi++) {
      double phi = phiMin + iphi*phiStep;

      n_dir.setRThetaPhi(1., theta, phi);		// Unit vector

      for (int mode = 0; mode < G4PhononPolarization::NUM_MODES; mode++) {
	vphase = mapper->getPhaseSpeed(mode, n_dir);
	vgroup = mapper->getGroupVelocity(mode, n_dir);
	slowness = mapper->getSlowness(mode, n_dir);
	polarization = mapper->getPolarization(mode, n_dir);

        lookupData[mode][N_X].push_back(n_dir.x());	// Wavevector dir.
	lookupData[mode][N_Y].push_back(n_dir.y());
	lookupData[mode][N_Z].push_back(n_dir.z());
	lookupData[mode][THETA].push_back(theta);	// ... and angles
	lookupData[mode][PHI].push_back(phi);
	lookupData[mode][S_X].push_back(slowness.x());	// Slowness direction
	lookupData[mode][S_Y].push_back(slowness.y());
	lookupData[mode][S_Z].push_back(slowness.z());
	lookupData[mode][S_MAG].push_back(slowness.mag());
	lookupData[mode][S_PAR].push_back(slowness.perp());
	lookupData[mode][V_P].push_back(vphase);
	lookupData[mode][V_G].push_back(vgroup.mag());
	lookupData[mode][V_GX].push_back(vgroup.x());
	lookupData[mode][V_GY].push_back(vgroup.y());
	lookupData[mode][V_GZ].push_back(vgroup.z());
	lookupData[mode][E_X].push_back(polarization.x());
	lookupData[mode][E_Y].push_back(polarization.y());
	lookupData[mode][E_Z].push_back(polarization.z());
      }
    }
  }
}

// *****************************************************************************

// $$$$$$$$$$$$$$$$$$$$$$$$$$$ EVEN INTERPOLATION HEADERS $$$$$$$$$$$$$$$$$$$$$$

/* given the (pointer to the) evenly spaced interpolation grid
   generated previously, this method returns an interpolated value for
   the data type already built into the G4CMPGridInterp structure */
double G4CMPPhononKinTable::interpolateEven(G4CMPGridInterp& grid,
					    double theta, double phi) {
    // check that the n values we're interpolating at are possible:
  if (!goodBin(theta,phi)) {
    cerr << "ERROR: Cannot interpolate (" << theta << ", " << phi << ")"
	 << endl;
    return ERRONEOUS_INPUT; // exit method
  }

  // perform interpolation and return result:
  return grid.interp(theta, phi);
}

/* overloaded method interpolates on the vector of vectors of
   interpolation grid structures.  Much the same as its simpler version,
   but this one requires specification of the mode and data type
   desired */
double G4CMPPhononKinTable::interpolateEven(double theta, double phi, int MODE,
					    int TYPE_OUT, bool SILENT) {
  // check that the n values we're interpolating at are possible:
  if (!goodBin(theta,phi)) {
    cerr << "ERROR: Cannot interpolate (" << theta << ", " << phi << ")"
	 << endl;
    return ERRONEOUS_INPUT; // exit method
  }
  
  // check that the data type we're asking for is not theta or phi
  if (TYPE_OUT == THETA || TYPE_OUT == PHI) {
    if (!SILENT) cout << "Warning: Data type desired provided as input\n";
    // no sense in actually interpolating here, just return appropriate input
    return (TYPE_OUT==THETA ? theta : phi);
  }
    
  // perform interpolation and return result:
  return interpolateEven(quantityMap[MODE][TYPE_OUT], theta, phi);
}

/* sets up JUST ONE interpolation table for N_X and N_Y, which are evenly spaced.
   Aphi one kind of data (TYPE_OUT) can be read off of this table */
G4CMPGridInterp 
G4CMPPhononKinTable::generateEvenTable(int MODE,
				     G4CMPPhononKinTable::DataTypes TYPE_OUT) {
  /* set up the two vectors of input components (N_X and N_Y) which
     will define the grid that will be interpolated on
     these grids are deliberately oversized, as our grid is evenly
     spaced but circular, not square. */
  vector<double> x1;
  for (int ith=0; ith<=thetaCount; ith++) x1.push_back(thetaMin+ith*thetaStep);

  vector<double> x2;
  for (int iphi=0; iphi<=phiCount; iphi++) x2.push_back(phiMin+iphi*phiStep);

  /* set up the matrix of data values to be interpolated */
  matrix<double> dataVals(x1.size(), x2.size(), OUT_OF_BOUNDS);
  size_t SIZE = lookupData[MODE][N_X].size();
  for (size_t i = 0; i<SIZE; i++) {
    double thetaVal = lookupData[MODE][THETA][i];
    double phiVal = lookupData[MODE][PHI][i];
    
    // figure out n_x-index of current data point:
    size_t x1vecIndex = (size_t)((thetaVal-thetaMin) / thetaStep);
    if (doubApproxEquals(x1[x1vecIndex+1], thetaVal) && x1vecIndex+1 < x1.size())
      x1vecIndex++;                           // resolve possible rounding issue
    else if (!doubApproxEquals(x1[x1vecIndex], thetaVal)) { // check
      cerr << "ERROR: Unanticipated index rounding behavior, thetaVal= "
	   << thetaVal << " index " << x1vecIndex << " = " << x1[x1vecIndex]
	   << endl;
    }

    // figure out n_y-index of current data point:
    size_t x2vecIndex = (size_t)((phiVal-phiMin) / phiStep);    // may round down to nearest integer index
    if (doubApproxEquals(x2[x2vecIndex+1], phiVal) && x2vecIndex+1 < x2.size())
      x2vecIndex++;                           // resolve possible rounding issue
    else if (!doubApproxEquals(x2[x2vecIndex], phiVal)) {  // check
      cerr << "ERROR: Unanticipated index rounding behavior, phiVal= "
	   << phiVal << " index " << x2vecIndex << " = " << x2[x2vecIndex]
	   << endl;
    }

    // fill in data point in matrix:
    // first check to make sure that point not already entered:
    if (!doubApproxEquals(dataVals[x1vecIndex][x2vecIndex], OUT_OF_BOUNDS,
			  1.0e-10, false))
      cerr << "ERROR: Attempting to enter a data value more than once @ "
	   << x1vecIndex << " " << x2vecIndex << endl;
    else // at long last, put data point in matrix
      dataVals[x1vecIndex][x2vecIndex] = lookupData[MODE][TYPE_OUT][i];
  }

  // create and return interpolation data structure, the output of this method
  return G4CMPGridInterp(x1, x2, dataVals);
}

/* a method for setting up a vector of interpolation grids - it does
   so by calling the single grid constructor maphi times. While this is
   certainly not the absolute most efficient method concievable, the
   gains in efficiency that might be made are not terribly
   significant. Output is a vector of vectors of G4CMPGridInterps, with the
   first index specifying mode, and the 2nd data type, as ususal */
void G4CMPPhononKinTable::generateMultiEvenTable() {
  clearQuantityMap();

  for (int mode = 0; mode < G4PhononPolarization::NUM_MODES; mode++) {
    // create 2nd level vectors:
    vector<G4CMPGridInterp> subTable;
    for (int dType = 0; dType < NUM_DATA_TYPES; dType++)
      // make individual tables:
      subTable.push_back(generateEvenTable(mode, (DataTypes)dType));
    quantityMap.push_back(subTable);
  }
}

// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

// +++++++++++++++++++++++++++++ COMPLETE LOOKUP TABLE +++++++++++++++++++++++++
void G4CMPPhononKinTable::write() {
  // <^><^><^><^><^><^><^><^><^> INITIAL SETUP <^><^><^><^><^><^><^><^><
  // set up the lookup table as a data file:
  string fName = mapper->getLatticeName()+"LookupTable.txt";
  ofstream lookupTable(fName.c_str());
  if (!lookupTable.good()) {
    G4ExceptionDescription msg;
    msg << "Unable to open " << fName << " for output.";
    G4Exception("G4CMPPhononKinTable::write", "Phonon010",
		FatalException, msg);
    return;
  }
  lookupTable.setf(ios::left);

  // file header:
  /* the next few steps ensure the same number of header lines looked for when reading
     from the lookup table */
  vector<string> headerLines;
  headerLines.push_back("Columns:");
  headerLines.push_back("mode");
  headerLines.push_back("----");
  for (int i = 0; i < NUM_DATA_TYPES; i++) {
    headerLines[1] += " | " + getDataTypeName(i);
    headerLines[2] += "------";
  }

  for (size_t i=0; i < headerLines.size(); i++)
    lookupTable << "# " << headerLines[i] << endl;
  // <^><^><^><^><^><^><^><^><^><^><^><^><^><^><^><^><^><^><^><^><^><^><

  int entry=0;
  for (int ith=0; ith<=thetaCount; ith++) {
    for (int iphi=0; iphi<=phiCount; iphi++) {
      for (int mode = 0; mode < G4PhononPolarization::NUM_MODES; mode++) {
	lookupTable << setw(18) << G4PhononPolarization::Label(mode);
	for (size_t cols=0; cols<NUM_DATA_TYPES; cols++) {
	  lookupTable << setw(18) << lookupData[mode][cols][entry]/getDataUnit(cols);
	}
	lookupTable << endl;
      }

      entry++;		// Increment counter, corresponding to push_back()
    }
  }
}

// given the data type index, returns the abbreviation (s_x, etc...)
// make sure to update this method if the data types are altered
string G4CMPPhononKinTable::getDataTypeName(int TYPE) {
  switch (TYPE) {
  case N_X:   return "n_x";
  case N_Y:   return "n_y";
  case N_Z:   return "n_z";
  case THETA: return "theta";
  case PHI:   return "phi";
  case S_X:   return "s_x";
  case S_Y:   return "s_y";
  case S_Z:   return "s_z";
  case S_MAG: return "s";
  case S_PAR: return "s_par";
  case V_P:   return "v_p";
  case V_G:   return "V_g";
  case V_GX:  return "V_gx";
  case V_GY:  return "V_gy";
  case V_GZ:  return "V_gz";
  case E_X:   return "e_x";
  case E_Y:   return "e_y";
  case E_Z:   return "e_z";
  default:  throw("ERROR: not al valid data type");
  } throw("ERROR: not al valid data type");
}

// given the data type index, returns the G4 units for writing output
// make sure to update this method if the data types are altered
G4double G4CMPPhononKinTable::getDataUnit(int TYPE) {
  switch (TYPE) {
  case N_X:   return 1.;	// Unit vector
  case N_Y:   return 1.;
  case N_Z:   return 1.;
  case THETA: return rad;
  case PHI:   return rad;
  case S_X:   return s/m;	// Inverse velocity
  case S_Y:   return s/m;
  case S_Z:   return s/m;
  case S_MAG: return s/m;
  case S_PAR: return s/m;
  case V_P:   return m/s;
  case V_G:   return m/s;
  case V_GX:  return m/s;
  case V_GY:  return m/s;
  case V_GZ:  return m/s;
  case E_X:   return 1.;	// Unit vector
  case E_Y:   return 1.;
  case E_Z:   return 1.;
  default:  throw("ERROR: not al valid data type");
  } throw("ERROR: not al valid data type");
}
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// ############################# EXTRANEOUS HEADERS ############################

/* takes two Doubs and returns 'true' if they are sufficiently close
   to one another.  This is to get around the often annoying issue of
   doubles not quite being equal, even though they should be, in an
   exact sense. Set TOL to 0 to do a normal '==' check */
bool doubApproxEquals(double D1, double D2, double TOL, bool USE_ABS_TOL)
{
    // determine value of tolerance given type of tolerance to use
    double absTol = TOL;                                  // use absolute tolerance
    if (!USE_ABS_TOL) absTol = ((D1 + D2)/2.0) * TOL;   // use relative tolerance
    
    // check for (near) equality:
    return (D1-D2 >= -absTol && D1-D2 <= absTol);
}

/* takes two Doubs and returns true if the first is significantly
   greater than the second, with "significantly" here meaning by an
   amount more than the tolerance. Set TOL to 0 to do a normal '>'
   check */
bool doubGreaterThanApprox(double D1, double D2, double TOL, bool USE_ABS_TOL)
{
    // determine value of tolerance given type of tolerance to use
    double absTol = TOL;                                  // use absolute tolerance
    if (!USE_ABS_TOL) absTol = ((D1 + D2)/2.0) * TOL;   // use relative tolerance
    
    // check for (near) equality:
    return (D1 > D2 + absTol);	                  // DEFINITELY greater than
}

/* takes two Doubs and returns true if the first is significantly less
   than the second, with "significantly" here meaning by an amount more
   than the tolerance Set TOL to 0 to do a normal '<' check*/
bool doubLessThanApprox(double D1, double D2, double TOL, bool USE_ABS_TOL)
{
    // determine value of tolerance given type of tolerance to use
    double absTol = TOL;                                  // use absolute tolerance
    if (!USE_ABS_TOL) absTol = ((D1 + D2)/2.0) * TOL;   // use relative tolerance
    
    // check for (near) equality:
    return (D1 < D2 - absTol);                  // DEFINITELY less than
}
// #############################################################################

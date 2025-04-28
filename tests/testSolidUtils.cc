// This script tests a full functionality of the transforms done in the
// SolidUtils class. This script will test all 4 ways to set the different
// transforms: null/Identity, Rotation and Translations,
// G4AffineTransform object, and a touchable.
//
// Usage: testSolidUtils x y z dx dy dz verboseLevel
//
// 20250220  M. Kelsey -- Create iZIP5 construction
// 20250428  N. Tenpas -- Create test for SolidUtils class.

#include "G4CMPSolidUtils.hh"
#include "G4AffineTransform.hh"
#include "G4Box.hh"
#include "G4IntersectionSolid.hh"
#include "G4PhysicalConstants.hh"
#include "G4RunManager.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#include "G4VSolid.hh"
#include "G4VTouchable.hh"

/*#include <float.h>
#include <cstdlib> 
#include <string>*/

G4VSolid* MakeZIP() {
  G4double ZipRad = 38.1*mm;		// iZIP5, 3" x 1"
  G4double ZipThick = 25.4*mm - 1.2*um;	// This makes space for QETs
  G4double ZipAxis1Len = 75.4888*mm;	// X-axis flats
  G4double ZipAxis2Len = 72.1868*mm;	// Y-axis flats

  G4ThreeVector ExtraFlats;		// This is the diagonal flat
  ExtraFlats.setRhoPhiZ(37.7444*mm, 45.*deg, 0.);

  G4double tolerance = 1.*nm;		// Little bit extra for subtractions

  G4String GetName = "iZIP5";

  // Start with a cylindrical crystal
  G4VSolid* zipShape = new G4Tubs(GetName, 0,ZipRad,ZipThick/2.,0,360*deg);

  // Make both side cuts in one step
  G4Box* zipCutBox = new G4Box("ZipCutBox", ZipAxis1Len/2., ZipAxis2Len/2.,
			       (ZipThick+tolerance)/2.);
  zipShape = new G4IntersectionSolid(GetName, zipShape, zipCutBox);

  // Remove additional single flats if necessary
  G4double phi   = ExtraFlats.phi();
  G4double rflat = ExtraFlats.r();
  if (rflat < ZipRad) {
    G4double delta = ZipRad - rflat;          	// Sagitta to be cut away
    G4double hwid  = sqrt(ZipRad*ZipRad - rflat*rflat);

    // Make box with double the sagitta wide, and double the chord long
    G4Box* flatCut = new G4Box("FlatCut", delta, hwid, (ZipThick+tolerance)/2.);

    G4ThreeVector pos;
    pos.setRhoPhiZ(ZipRad, phi, 0.);		// Place on rim of cylinder

    zipShape = new G4SubtractionSolid(GetName, zipShape, flatCut,
				      new G4RotationMatrix(phi,0.,0.), pos);
  }

  return zipShape;
}


int main(int argc, char** argv) {
  if (argc < 8) { return 1; }
  G4double X = std::atof(argv[1]);
  G4double Y = std::atof(argv[2]);
  G4double Z = std::atof(argv[3]);
  G4ThreeVector pos(X, Y, Z);

  G4double dx = std::atof(argv[4]);
  G4double dy = std::atof(argv[5]);
  G4double dz = std::atof(argv[6]);
  G4ThreeVector dir(kx, ky, kz);

  G4int verboseLevel = std::atoi(argv[7]);

  G4RunManager* runManager = new G4RunManager();
  G4VSolid* solid = MakeZIP();

  // Test 1: Position and direction transforms with identity
  G4CMPSolidUtils* solidUtils = new G4CMPSolidUtils(solid, verboseLevel,
                                                    "Test");
  G4cout << "TEST 1:\n  Original Position: " << pos
    << "\n  Original Direction: " << dir << G4endl;

  // Test 3: Position transforms with given transform

  // Test 4: Direction transforms with given transform

  // Test 5: Position transforms with given G4AffineTransform

  // Test 6: Direction transforms with given G4AffineTransform

  // Test 7: Position transform with touchable

  // Test 8: Direction transform with touchable

  delete runManager;
  return 0;
}

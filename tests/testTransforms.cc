// This script tests a full functionality of the transforms done in the
// SolidUtils class. This script will test 3 ways to set the different
// transforms: null/Identity, with rotation and translation matrices,
// and with a G4AffineTransform object.
//
// Usage: testSolidUtils x y z dx dy dz verboseLevel
//
// Where: The position to be transformed will be (x, y, z)
//        The direction to be transformed will be (dx, dy, dz)
//
// 20250220  M. Kelsey -- Create iZIP5 construction
// 20250505  N. Tenpas -- Create transform test

#include "G4AffineTransform.hh"
#include "G4Box.hh"
#include "G4IntersectionSolid.hh"
#include "G4PhysicalConstants.hh"
#include "G4RotationMatrix.hh"
#include "G4RunManager.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#include "G4VSolid.hh"
#include "G4VTouchable.hh"


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


int main(/*int argc, char* argv[]*/) {
  G4ThreeVector dir(1, 0, 0);
  G4ThreeVector origDir = dir;
  G4RunManager* runManager = new G4RunManager();

  // Rotation Matrix to transform myDir
  G4RotationMatrix rotM = G4RotationMatrix::IDENTITY;
  G4double rotAng = 45 * deg;
  rotM.rotateZ(rotAng);

  // G4AffineTransform Object
  G4ThreeVector transV(0,0,0);
  G4AffineTransform transObj = G4AffineTransform(rotM, transV);

  // -------------------------------------------------------
  G4cout << "Rotation Matrix: " << rotM << "\nAffine Object:\n " << transObj << "\n" << G4endl;

  G4cout << "Original Direction: " << dir << G4endl;
  
  // Lattice transform
  G4ThreeVector testDir = rotM * dir;
  G4ThreeVector testDir2 = rotM.inverse() * dir;
  G4cout << "New Direction with v.transform(ROT): " << testDir << G4endl;
  G4cout << "Inverse Direction with v.transform(ROT): " << testDir2 << "\n" << G4endl;

  // Geom Utils Transform
  G4ThreeVector affDir = transObj.TransformAxis(dir);
  G4cout << "Direction with TransformAxis: " << affDir << "\n" << G4endl;

  // Solit Utils Transform
  transObj.ApplyAxisTransform(dir);
  G4cout << "Direction with ApplyAxisTransform: " << dir << G4endl;

  delete runManager;
  return 0;
}

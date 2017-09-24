/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// Usage: testChargeCloud <N> <Lattice>
//
// Specify number of point to throw, and which lattice directory to use.
// Geant4 material will be set as "G4_<Lattice>".
//
// NOTE: 10 keV energy deposit should produce ~5000 e/h pairs

#include "globals.hh"
#include "G4CMPChargeCloud.hh"
#include "G4LatticeManager.hh"
#include "G4LatticeLogical.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include <stdlib.h>
#include <vector>


int main(int argc, char* argv[]) {
  if (argc < 3) {
    G4cerr << "Usage: " << argv[0] << " <N> <Lattice>" << G4endl;
    ::exit(1);
  }

  G4int nErrors = 0;			// Increment counter at sanity checks

  G4int npoints = atoi(argv[1]);
  G4String lname = argv[2];
  G4String mname = "G4_"+lname;

  G4Material* mat = G4NistManager::Instance()->FindOrBuildMaterial(mname);
  G4LatticeLogical* lat = G4LatticeManager::Instance()->LoadLattice(mat,lname);

  // MUST USE 'new', SO THAT G4SolidStore CAN DELETE
  G4Tubs* crystal = new G4Tubs("GeCrystal", 0., 5.*cm, 1.*cm, 0., 360.*deg);

  G4CMPChargeCloud cloud(lat, crystal);
  cloud.SetVerboseLevel(2);

  G4double rcloud = cloud.GetRadius(npoints);

  G4cout << "G4CMPChargeCloud " << npoints << " e/h radius " << rcloud/nm
	 << " nm" << G4endl;

  std::vector<G4ThreeVector> points;

  // Generate points in bulk of crystal (no boundary effects)
  G4ThreeVector center(0.,0.,0.);
  points = cloud.Generate(npoints, center);

  G4cout << " generated " << points.size() << " points" << G4endl;

  if (points.size() != (size_t)npoints) {
    G4cerr << " WRONG NUMBER OF POINTS" << G4endl;
    nErrors++;
  }

  // Find boundaries of good sphere for validation
  G4ThreeVector min(1.*km,1.*km,1.*km), max(-1.*km,-1.*km,-1*km);
  for (size_t i=0; i<points.size(); i++) {
    if (points[i].x() < min.x()) min.setX(points[i].x());
    if (points[i].y() < min.y()) min.setY(points[i].y());
    if (points[i].z() < min.z()) min.setZ(points[i].z());
    if (points[i].x() > max.x()) max.setX(points[i].x());
    if (points[i].y() > max.y()) max.setY(points[i].y());
    if (points[i].z() > max.z()) max.setZ(points[i].z());
  }

  G4cout << " at " << center << " points span " << min/nm << " to "
	 << max/nm << " nm" << G4endl;

  if (min.x() < -rcloud || min.y() < -rcloud || min.z() < -rcloud ||
      max.x() > rcloud || max.y() > rcloud || max.z() > rcloud) {
    G4cerr << " POINTS GENERATED OUTSIDE RADIUS" << G4endl;
    nErrors++;
  }

  // Generate points close to corner (top face and side effects)
  G4ThreeVector edge(0., 5.*cm - rcloud/2., 1.*cm - rcloud/2.);
  points = cloud.Generate(npoints, edge);

  // Find boundaries of good sphere for validation
  min.set(1.*km,1.*km,1.*km); max.set(-1.*km,-1.*km,-1*km);
  for (size_t i=0; i<points.size(); i++) {
    if (points[i].x() < min.x()) min.setX(points[i].x());
    if (points[i].y() < min.y()) min.setY(points[i].y());
    if (points[i].z() < min.z()) min.setZ(points[i].z());
    if (points[i].x() > max.x()) max.setX(points[i].x());
    if (points[i].y() > max.y()) max.setY(points[i].y());
    if (points[i].z() > max.z()) max.setZ(points[i].z());
  }

  G4cout << " at " << center << " points span " << min/nm << " to "
	 << max/nm << " nm" << G4endl;
  if (max.x() > 5.*cm || max.y() > 5.*cm || max.z() > 1.*cm) {
    G4cerr << " POINTS GENERATED OUTSIDE VOLUME" << G4endl;
    nErrors++;
  }

  ::exit(nErrors);
}

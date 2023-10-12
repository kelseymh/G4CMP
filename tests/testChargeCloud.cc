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

// Global variables for use in tests

namespace {
  G4CMPChargeCloud* cloud = 0;
  G4int nErrors = 0;		// Increment counter at failed checks
}


// Test cloud generation at specific location

void testCloud(G4int n, const G4ThreeVector& pos) {
  const std::vector<G4ThreeVector>& points = cloud->Generate(n, pos);

  G4cout << " generated " << points.size() << " points" << G4endl;

  if (points.size() != (size_t)n) {
    G4cerr << " WRONG NUMBER OF POINTS" << G4endl;
    nErrors++;
  }

  // Find boundaries of good sphere for validation
  G4ThreeVector min(1.*km,1.*km,1.*km);
  G4ThreeVector max(-1.*km,-1.*km,-1*km);
  G4int maxbin = -1;
  G4double rsum=0., r2sum=0.;

  for (size_t i=0; i<points.size(); i++) {
    G4int ibin = cloud->GetPositionBin(i);
    if (ibin > maxbin) maxbin = ibin;

    G4double ri = (points[i]-pos).mag();
    rsum += ri;
    r2sum += ri*ri;

    if (points[i].x() < min.x()) min.setX(points[i].x());
    if (points[i].y() < min.y()) min.setY(points[i].y());
    if (points[i].z() < min.z()) min.setZ(points[i].z());
    if (points[i].x() > max.x()) max.setX(points[i].x());
    if (points[i].y() > max.y()) max.setY(points[i].y());
    if (points[i].z() > max.z()) max.setZ(points[i].z());
  }

  G4double ravg = rsum/points.size();
  G4double rrms = sqrt(r2sum/points.size() - ravg*ravg);
  G4cout << " at " << pos << " points span " << min/mm << " to "
	 << max/mm << " mm\n Ravg " << ravg/nm << " rms " << rrms/nm << " nm"
	 << "\n Maximum bin (zzzyyyxxx) " << maxbin << G4endl;

  G4double rcloud = cloud->GetRadius();		// Radius used to generate

  if ((min-pos).x() < -rcloud ||
      (min-pos).y() < -rcloud ||
      (min-pos).z() < -rcloud ||
      (max-pos).x() > rcloud ||
      (max-pos).y() > rcloud ||
      (max-pos).z() > rcloud) {
    G4cerr << " POINTS GENERATED OUTSIDE RADIUS" << G4endl;
    nErrors++;
  }
}


// Main test is here

int main(int argc, char* argv[]) {
  if (argc < 3) {
    G4cerr << "Usage: " << argv[0] << " <N> <Lattice> [verbose]" << G4endl;
    ::exit(1);
  }

  G4int npoints = atoi(argv[1]);
  G4String lname = argv[2];
  G4String mname = "G4_"+lname;

  G4int verbose = (argc>3) ? atoi(argv[3]) : 0;

  G4Material* mat = G4NistManager::Instance()->FindOrBuildMaterial(mname);
  G4LatticeLogical* lat = G4LatticeManager::Instance()->LoadLattice(mat,lname);

  // MUST USE 'new', SO THAT G4SolidStore CAN DELETE
  G4Tubs* crystal = new G4Tubs("GeCrystal", 0., 5.*cm, 1.*cm, 0., 360.*deg);

  cloud = new G4CMPChargeCloud(lat, crystal);
  cloud->SetVerboseLevel(verbose);

  G4double rcloud = cloud->ComputeRadius(npoints);

  G4cout << "G4CMPChargeCloud " << npoints << " e/h radius "
	 << rcloud/nm << " nm" << G4endl;

  // Generate points in bulk of crystal (no boundary effects)
  testCloud(npoints, G4ThreeVector(0.,0.,0.));

  // Generate points close to corner (top face and side effects)
  testCloud(npoints, G4ThreeVector(0., 5.*cm-rcloud/2., 1.*cm-rcloud/2.));

  delete cloud;		// Clean up memory at end
  ::exit(nErrors);
}

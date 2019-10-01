/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/
/*
 * A test for checking changes to the EField code.
 *
 * Usage: g4cmpEFieldText lx ly lz N [3D|2D]
 *
 * Takes in the x, y, and z dimensions of a rectangular box with a number
 * of points N to sample in the box, and whether to test with a full 3D
 * or 2D projective (r-z) mesh.
 *
 * Outputs a file with the voltage and E field components at each position.
 * Also prints how long it takes to run.
 *
 * 20170527  Abort job if output file fails
 * 20190918  Convert to test either 2D or 3D mesh; input names predefined
 */

#include "globals.hh"
#include "geomdefs.hh"
#include "G4CMPMeshElectricField.hh"
#include "G4Tokenizer.hh"
#include <stdlib.h>
#include <array>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
using namespace std;


void testField(G4double lx, G4double ly, G4double lz, G4int N,
	       G4CMPMeshElectricField& field, const G4String& outname) {
  ofstream outputFile(outname);
  if (!outputFile.is_open()) {
    cerr << "Cannot create " << outname << ": " << strerror(errno)
	      << endl;
    ::exit(1);
  }

  outputFile << "   x\tt   y\tt   z\tt   V\t\t   Ex\t\t   Ey\t\t   Ez"
	     << endl;

  G4double EMfield[6];
  G4int n = cbrt(N);
  G4double Deltax = lx/n;
  G4double Deltay = ly/n;
  G4double Deltaz = lz/n;

  G4cout << "Generating " << n*n*n << " points in steps of " << Deltax
	 << " " << Deltay << " " << Deltaz << " ..." << G4endl;

  for (G4int i = 0; i < n; ++i) {
    for (G4int j = 0; j < n; ++j) {
      for (G4int k = 0; k < n; ++k) {
	G4double x = -(Deltax * n/2) + (i * Deltax);
	G4double y = -(Deltay * n/2) + (j * Deltay);
	G4double z = -(Deltaz * n/2) + (k * Deltaz);
	G4double pt[4] = {x, y, z, 0};
	G4double V = field.GetPotential(pt);
	field.GetFieldValue(pt, EMfield);

	outputFile << x << "\t" << y << "\t" << z << "\t" << V << "\t"
		   << EMfield[3] << "\t" << EMfield[4] << "\t" << EMfield[5]
		   << endl;
      }
    }
  }
  
  outputFile.close();
}


// Special function to find pure numeric-data lines

G4bool isNumbers(const G4String& line, size_t len=0) {
  size_t nonnum = line.find_first_not_of("0123456789.Ee+- \t");
  return (nonnum >= (len==0 ? line.length() : len));
}


// Note that input filenames are predefined

void readMesh(const G4String& type, vector<array<G4double,3> >& xyz,
	      vector<G4double>& voltage, vector<array<G4int,4> >& tetra) {
  xyz.clear();
  voltage.clear();
  tetra.clear();

  G4String line;	// Reusable input buffers
  array<G4double, 3> pos;
  array<G4int, 4> simplex;

  G4String pfile = "EField"+type+"_points.dat";

  ifstream points(pfile);
  if (!points.good()) {
    G4cerr << " Error reading " << pfile << G4endl;
    ::exit(2);
  }

  while (points.good()) {
    getline(points, line);
    if (line.empty() || line(0) == '#' || !isNumbers(line, 7)) continue;

    G4Tokenizer values(line);
    pos[0] = strtod(values().c_str(), NULL);
    pos[1] = strtod(values().c_str(), NULL);
    pos[2] = (type=="3D" ? strtod(values().c_str(), NULL) : 0.);
    xyz.push_back(pos);
    voltage.push_back(strtod(values().c_str(), NULL));
  }
  points.close();

  G4String tfile = "EField"+type+"_tetra.dat";
  ifstream tetras(tfile);
  if (!tetras.good()) {
    G4cerr << " Error reading " << tfile << G4endl;
    ::exit(2);
  }

  while (tetras.good()) {
    getline(tetras, line);
    if (line.empty() || line(0) == '#' || !isNumbers(line, 1)) continue;

    G4Tokenizer values(line);
    simplex[0] = atoi(values().c_str());
    simplex[1] = atoi(values().c_str());
    simplex[2] = atoi(values().c_str());
    simplex[3] = (type=="3D" ? atoi(values().c_str()) : -1);
    tetra.push_back(simplex);
  }
  tetras.close();
}


// Driver program for testing

int main(int argc, char* argv[]) {
  if (argc != 6) {
    G4cout << argv[0] << " xlen ylen zlen Npoints 3D|2D" << G4endl;
    return 1;
  }

  G4double lx = strtod(argv[1], NULL);		// Command line arguments
  G4double ly = strtod(argv[2], NULL);
  G4double lz = strtod(argv[3], NULL);
  G4int N = atoi(argv[4]);
  G4String type = argv[5];

  if (type != "3D" && type != "2D") {
    G4cerr << "Must specify 3D or 2D testing." << G4endl;
    return 1;
  }

  // Read predefined mesh files for 2D or 3D testing
  vector<array<G4double,3> > xyz;
  vector<G4double> voltage;
  vector<array<G4int,4> > tetra;

  readMesh(type, xyz, voltage, tetra);
  G4cout << type << " read " << xyz.size() << " points, "
	 << tetra.size() << " simplices" << G4endl;

  // 2D projection is (r-z), 3D is full coordinates
  EAxis xproj = (type=="2D" ? kRho : kUndefined);
  EAxis yproj = (type=="2D" ? kZAxis : kUndefined);

  time_t start, end;

  time(&start);
  G4CMPMeshElectricField EField(xyz, voltage, tetra, xproj, yproj);
  testField(lx, ly, lz, N, EField, "EFieldTest.txt");
  time(&end);
  G4cout << type << " test took " << difftime(end, start) << " s" << G4endl;

  return 0;
}

/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/
/*
 * A test for checking changes to the EField code.
 * Takes in the x, y, and z dimensions of a rectangular box, the
 * number of points, and an EPot file name. Uniformly distributes
 * points in the rectangular box centered at the origin using the
 * given potential from the EPot file.
 * Outputs a file with the voltage and E field components at each position.
 * Also prints how long it takes to run.
 *
 * 20170527  Abort job if output file fails
 * 20180712  Expand to exercise field manager, different field types
 */

#include "G4CMPMeshElectricField.hh"
#include "G4CMPFieldManager.hh"
#include "G4CMPFieldUtils.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4UniformElectricField.hh"
#include <ctime>
#include <fstream>
#include <stdlib.h>


// Validate mesh file input

G4CMPMeshElectricField* 
testMeshFile(const G4String& filename, G4int nstep,
	     G4double Deltax, G4double Deltay, G4double Deltaz) {
  G4CMPMeshElectricField* EField = new G4CMPMeshElectricField(filename);

  std::fstream outputFile("EFieldTestdata.txt", std::fstream::out);
  if (!outputFile.is_open()) {
    G4cerr << "Cannot create file, error was: " << strerror(errno)
	  << G4endl;
    return EField;
  }
  outputFile << "Data formatted as x, y, z, V, Ex, Ey, Ez"
	     << std::endl;

  time_t start, end;
  std::time(&start);

  G4double EB[6] = {0,0,0,0,0,0};	// Field buffer: Bx,By,Bz,Ex,Ey,Ez
  G4double pos[4] = {0,0,0,0};		// Position buffer: x,y,z,t

  for (G4int i = 0; i < nstep; ++i) {
    pos[0] = i*Deltax - Deltax*nstep/2;

    for (G4int j = 0; j < nstep; ++j) {
      pos[1] = j*Deltay - Deltay*nstep/2;

      for (G4int k = 0; k < nstep; ++k) {
	pos[2] = k*Deltaz - Deltaz*nstep/2;
	EField->GetFieldValue(pos, EB);

	outputFile << pos[0] << " " << pos[1] << " " << pos[2] << " "
		   << EField->GetPotential(pos) << " "
		   << EB[3] << " " << EB[4] << " " << EB[5] << " "
		   << std::endl;
      }
    }
  }

  outputFile.close();

  std::time(&end);
  G4cout << "testMeshFile took " << difftime(end, start) << " s" << G4endl;
  
  return EField;
}


// Exercise G4CMPFieldUtils at fixed position in mesh

void testFieldUtils(G4CMPMeshElectricField* field, G4double x, G4double y,
		    G4double z) {
  G4Tubs* crystalS = new G4Tubs("Crystal", 0., 75.*mm, 12.7*mm, 0., 360.*deg);
  G4Material* Ge = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ge");
  G4LogicalVolume* crystal = new G4LogicalVolume(crystalS, Ge, "Crystal");
  crystal->SetFieldManager(new G4CMPFieldManager(field), false);

  G4ThreeVector pos(x,y,z);

  /*** DON'T HAVE (vol,pos) UTILITIES ANY MORE.  WHAT TO DO?
  G4cout << "testFieldUtils interrogating @ " << pos
	 << "\n field vector " << G4CMP::GetFieldAtPosition(crystal,pos)/(volt/m)
	 << " = " << G4CMP::GetFieldAtPosition(crystal,pos).mag()/(volt/m)
	 << " V/m" << G4endl;

  G4cout << " potential " << G4CMP::GetPotentialAtPosition(crystal,pos)/volt
	 << " volt " << G4endl;
  ***/
}


// Main driver program

int main(int argc, char** argv) {
  if (argc != 6) {
    G4cout << argv[0] << " x_length y_length z_length number_points EPotfile"
	   << G4endl;
    return 0;
  }

  G4double lx=strtod(argv[1],NULL);
  G4double ly=strtod(argv[2],NULL);
  G4double lz=strtod(argv[3],NULL);
  G4int N = atoi(argv[4]);
  G4String filename = argv[5];

  // Load mesh file and test interpolation in grid
  G4int n = cbrt(N);
  G4double dx = lx/n, dy = ly/n, dz = lz/n;
  G4CMPMeshElectricField* EField = testMeshFile(filename, n, dx, dy, dz);

  // Attach field to volume and test accessor function
  testFieldUtils(EField, dx*n/6., -dy*n/6., dz*n/3.);

  return 0;
}

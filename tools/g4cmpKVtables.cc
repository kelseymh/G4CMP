//
//  g4cmpKVtables -- Generate wavevector mapping tables for G4CMP
//
//  Uses G4CMPPhononKVgMap, written by Daniel Palken 2014.
//

#include "G4CMPPhononKVgMap.hh"
#include "G4ThreeVector.hh"
#include <fstream>
#include <iomanip>
#include <string>
#include <math.h>
using namespace std;


// Units of density are kg/m3, of Cij are Pa
// NOTE:  To reproduce D. Brandt's Vg tables, use Ge 5196.2792 kg/m3
Material Ge("Ge", 5323.0, 1.26e11, 0.44e11, 0.67e11);
Material Si("Si", 2329.0, 1.656e11, 0.639e11, 0.795e11);

int main(int argc, const char * argv[])
{
  // 1. choose material
  Material *materialToUse = 0;
  if (argc > 1) {
    if (string(argv[1]) == "Ge") materialToUse = &Ge;
    if (string(argv[1]) == "Si") materialToUse = &Si;
  } else materialToUse = &Ge;		// No command line arguments

  if (!materialToUse) {
    cerr << argv[0] << " Invalid material " << argv[1] << " specified" << endl;
    ::exit(1);
  }

  // 2. set up the entire interpolation infrastrcuture
  G4CMPPhononKVgMap *map = new G4CMPPhononKVgMap(materialToUse);
  
  map->writeLookupTable();		// Dump Dan's version of data file
  
  // 3. G4CMP uses (theta,phi) binning
  const int Ntheta = 161;		// These must match config.txt values
  const int Nphi = 321;
  
  // 4. Loop over modes (L, ST, FT) to set up output files
  for (int mode=G4CMPPhononKVgMap::L; mode<G4CMPPhononKVgMap::NUM_MODES; mode++) {
    string vgname = map->getModeName(mode) + ".ssv";
    ofstream vgfile(vgname, ios::trunc);
    vgfile << scientific << setprecision(7);
    
    string vdirname = map->getModeName(mode) + "Vec.ssv";
    ofstream vdirfile(vdirname, ios::trunc);
    vdirfile << scientific << setprecision(7);
    
    cout << "Generating " << materialToUse->getName() << " "
	 << map->getModeName(mode) << " files " << Ntheta << " x " << Nphi
	 << endl;
    
    G4ThreeVector kvec, Vg;
    double theta, phi;
    for (int itheta = 0; itheta<Ntheta; itheta++) {
      cout << "." << flush;

      theta = M_PI*itheta / Ntheta;		// Bin edges

      for (int iphi = 0; iphi<Nphi; iphi++) {
	phi = 2.*M_PI*iphi / Nphi;		// Bin edges

	// Convert phi, theta bin center to unit vector
	kvec.setRThetaPhi(1.,theta,phi);
	
	// 5. Get group velocity for given direction, write to file
	Vg = map->getGroupVelocity(mode, kvec);
	
	vgfile << setw(16) << Vg.mag() << endl;
	
	vdirfile << setw(16) << Vg.x() << setw(16) << Vg.y()
		 << setw(16) << Vg.z() << endl;
      }	// phi
    }	// theta
    cout << endl;
  }		// mode
}

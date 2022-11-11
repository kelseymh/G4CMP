/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/* This test creates slowness surface and group velocity data for silicon.
 * Use plot_test_phonon_kinematics.py to plot the results.
 *
 * 20170527  Abort job if output files can't be opened
 * 20170620  Change 'is_good()' to 'good()'
 * 20221102  Fix units of slowness output; expand to write phase velocities.
 *		Generate points in cos(theta) steps, not theta.
 */

#include "G4CMPPhononKinematics.hh"
#include "G4LatticeLogical.hh"
#include "G4LatticeManager.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PhononPolarization.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <fstream>
#include <assert.h>

using CLHEP::pi;

void print_usage() {
    G4cout << "Usage: phononKinematics <Ge|Si>" << G4endl;
}

void useG4CMPSolver(G4LatticeLogical* lattice) {
  // Configure arrays indexed on phonon mode to make
  G4String suffix[] = { "long", "trans_slow", "trans_fast" };
  G4int nmode = G4PhononPolarization::NUM_MODES;

  G4String veltype[] = { "group_vel", "phase_vel", "slowness" };
  G4double velunit[] = { m/s,         m/s,         s/m        };
  G4int nvel = sizeof(veltype)/sizeof(G4String);

  std::ofstream veldata;	// Reuse this for each output file

  const G4double nbin = 100.;	// Number of angular bins to store
  
  G4CMPPhononKinematics solver(lattice);
  G4ThreeVector kdir, vel;

  // Loop over modes and velocity types (entries in switch MUST MATCH veltype)
  for (G4int mode=0; mode<nmode; mode++) {
    for (G4int vtype=0; vtype<nvel; vtype++) {
      G4String fname = lattice->GetName()+"_phonon_"+veltype[vtype]+"_"+suffix[mode];
      veldata.open(fname);
      assert(veldata.good());

      G4cout << "Filling " << fname << " ..." << G4endl;

      // Uniform steps in theta, adjust phi steps for uniform display
      for (G4double theta=0.; theta < halfpi; theta += halfpi/nbin) {
	G4double phibin = ceil(nbin*sin(theta));	// More phi at bottom
	for (G4double phi=0; phi < halfpi; phi += halfpi/phibin) {
	  if (theta == 0. && phi > 0.) break;	// Only need one point at pole
	  
	  kdir.setRThetaPhi(1., theta, phi);

	  switch (vtype) {	// These MUST MATCH veltype name order
	  case 0: vel = solver.getGroupVelocity(mode, kdir); break;
	  case 1: vel = solver.getPhaseSpeed(mode, kdir)*kdir; break;
	  case 2: vel = solver.getSlowness(mode, kdir); break;
	  default: G4cerr << "Invalid vtype " << vtype << G4endl; ::exit(1);
	  }

	  vel /= velunit[vtype];	// Prepare with units for output file

	  veldata << vel.x() << ", " << vel.y() << ", " << vel.z() << G4endl;
	}	// for (G4double phi
      }		// for (G4double theta

      veldata.close();
    }		// for (G4int vtype
  }		// for (G4int mode
}

int main(int argc, char** argv) {
  if (argc != 2) {
    print_usage();
    return 0;
  }

  const G4String matname = argv[1];
  G4Material* mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_"+matname);
  assert(mat);

  G4LatticeLogical* lattice = G4LatticeManager::Instance()->LoadLattice(mat, matname);
  assert(lattice);

  useG4CMPSolver(lattice);

  delete lattice;
  return 0;
}

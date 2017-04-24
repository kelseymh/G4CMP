/* This test creates slowness surface and group velocity data for silicon.
 * Use plot_test_phonon_kinematics.py to plot the results.
 */

#include "G4CMPPhononKinematics.hh"
#include "G4LatticeLogical.hh"
#include "G4LatticeReader.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4SystemOfUnits.hh"
#include <fstream>

using CLHEP::pi;

void print_usage() {
    G4cout << "Usage: phononKinematics <path to Si/config.txt>" << G4endl;
}

void useG4CMPSolver(G4LatticeLogical* lattice) {
  G4CMPPhononKinematics solver(lattice);

  std::ofstream slowness_trans_slow("phonon_slowness_trans_slow");
  std::ofstream slowness_trans_fast("phonon_slowness_trans_fast");
  std::ofstream slowness_longi("phonon_slowness_longi");
  std::ofstream trans_slow("phonon_group_vel_trans_slow");
  std::ofstream trans_fast("phonon_group_vel_trans_fast");
  std::ofstream longi("phonon_group_vel_long");

  G4ThreeVector kdir(1., 0, 0);
  for (G4double theta = 0; theta < pi / 2.; theta += pi / 200.) {
    for (G4double phi = 0; phi < pi / 2.; phi += pi / 200.) {
      kdir.setRThetaPhi(1., theta, phi);

      G4ThreeVector vg =
        solver.getGroupVelocity(G4PhononPolarization::TransFast, kdir);
      G4ThreeVector vp =
        solver.getSlowness(G4PhononPolarization::TransFast, kdir);
      trans_fast << vg.x()/m*s << ", "
                 << vg.y()/m*s << ", "
                 << vg.z()/m*s << G4endl;
      slowness_trans_fast << vp.x()/m*s << ", "
                          << vp.y()/m*s << ", "
                          << vp.z()/m*s << G4endl;

      vg = solver.getGroupVelocity(G4PhononPolarization::TransSlow, kdir);
      vp = solver.getSlowness(G4PhononPolarization::TransSlow, kdir);
      trans_slow << vg.x()/m*s << ", "
                 << vg.y()/m*s << ", "
                 << vg.z()/m*s << G4endl;
      slowness_trans_slow << vp.x()/m*s << ", "
                          << vp.y()/m*s << ", "
                          << vp.z()/m*s << G4endl;

      vg = solver.getGroupVelocity(G4PhononPolarization::Long, kdir);
      vp = solver.getSlowness(G4PhononPolarization::Long, kdir);
      longi << vg.x()/m*s << ", "
            << vg.y()/m*s << ", "
            << vg.z()/m*s << G4endl;
      slowness_longi << vp.x()/m*s << ", "
                     << vp.y()/m*s << ", "
                     << vp.z()/m*s << G4endl;
    }
  }

  slowness_trans_fast.close();
  slowness_trans_slow.close();
  slowness_longi.close();
  trans_fast.close();
  trans_slow.close();
  longi.close();
}

int main(int argc, char** argv) {
  if (argc != 2) {
    print_usage();
    return 0;
  }

  const G4String filename = argv[1];
  G4Material* silicon = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");

  G4LatticeLogical* lattice = G4LatticeReader().MakeLattice(filename);
  lattice->SetDensity(silicon->GetDensity());
  lattice->Initialize();

  useG4CMPSolver(lattice);

  delete lattice;
  return 0;
}

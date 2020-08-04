/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// electron_Epv.cc	Generate data showing electron kinematics, both
//			true valley-frame and G4 effective-mass.
//
// 20141216  Michael Kelsey
// 20200604  G4CMP-208: Eliminate "unused variable" warnings with output flag

#include "G4LatticeLogical.hh"
#include "G4LatticeManager.hh"
#include "G4Material.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4VelocityTable.hh"
#include <iostream>
#include <iomanip>
using namespace std;

const G4bool longOutput=false;


int main(/*int argc, char* argv[]*/) {
  G4Material* ge = new G4Material("Ge", 32., 72.630*g/mole, 5.323*g/cm3,
				  kStateSolid);
  G4LatticeLogical* lattice =
    G4LatticeManager::GetLatticeManager()->LoadLattice(ge, "Ge");

  // Electron wavevector corresponding to phonon speed in lattice
  G4double me_HV = lattice->GetElectronMass();
  G4double ksound = lattice->GetSoundSpeed() * me_HV / hbar_Planck;

  // Lookup table used by G4Track to compute velocity
  G4VelocityTable* velTable = G4VelocityTable::GetVelocityTable();

  // Buffers for computations inside loops
  G4double kel, Ekin, meff, vtrue, veff, vG4Track, EoverM;
  G4ThreeVector k_elec, p_elec, v_elec, kdir, pdir, vdir;

  // Set up output file for import to Excel
  if (longOutput) {
    cout << "kmag\t\tk/ks\tp(eV)\tth(deg)\t E(ueV)\tmeff\tv(km/s)\tveff\tvG4"
	 << endl;
  } else {
    cout << "Kx\tKy\tKz\tPx\tPy\tPz\tVx\tVy\tVz" << endl;
  }

  // Loop over wavevectors up to "Mach 4"
  for (G4double mach=0.05; mach<=4.; mach += 0.05) {
    kel = mach*ksound;

    // Loop over angles relative to valley
    for (G4double thdeg=0; thdeg<=90.; thdeg += 5.) {
      for (G4double phideg=0; phideg<360.; phideg += 5.) {
	if (phideg > 0. && thdeg==0.) break;

	k_elec.setRThetaPhi(kel, thdeg*deg, phideg*deg);
	p_elec = lattice->MapK_valleyToP(1, k_elec);
	v_elec = lattice->MapPtoV_el(1, p_elec);

	// Use momentum and valley to get energy, effective mass
	Ekin = lattice->MapPtoEkin(1, p_elec);
	meff = lattice->GetElectronEffectiveMass(1, p_elec);
	
	vtrue = v_elec.mag();			// True velocity
	veff = p_elec.mag() / meff / c_light;	// Non-relativistic
	
	// Compute velocity by copying G4Track method
	EoverM = Ekin/meff/c_squared;
	if (EoverM > velTable->GetMaxTOfVelocityTable()) {
	  vG4Track = c_light;
	} else if (EoverM < DBL_MIN) {
	  vG4Track =0.;
	} else if (EoverM < velTable->GetMinTOfVelocityTable()) {
	  vG4Track = c_light*std::sqrt(EoverM*(EoverM+2.))/(EoverM+1.0);
	} else {    
	  vG4Track = velTable->Value(EoverM);
	}
	
	// Report kinematics
	if (longOutput) {
	  cout << fixed
	       << setw(8) << setprecision(2) << kel << "\t" << mach
	       << "\t" << setprecision(4) << p_elec.mag()/(1e-6*MeV)
	       << "\t" << setw(2) << setprecision(0) << thdeg
	       << "\t" << setw(7) << setprecision(3) << Ekin*1e12/MeV
	       << "\t" << setprecision(5) << meff*c_squared/electron_mass_c2
	       << "\t" << setprecision(4) << vtrue/(km/s)
	       << "\t" << setprecision(4) << veff/(km/s)
	       << "\t" << setprecision(4) << vG4Track/(km/s) << endl;
	} else {
	  kdir = k_elec.unit();
	  pdir = p_elec.unit();
	  vdir = v_elec.unit();
	  cout << fixed << setprecision(4)
	       << kdir.x() << "\t" << kdir.y() << "\t" << kdir.z()
	       << "\t" << pdir.x() << "\t" << pdir.y() << "\t" << pdir.z()
	       << "\t" << vdir.x() << "\t" << vdir.y() << "\t" << vdir.z()
	       << endl;
	}
      }	/* phideg */
    }	/* thdeg */
  }	/* mach */
}

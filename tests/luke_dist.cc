/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// luke_dist.cc		Generate data showing how phonon HV vector varies
//			with k and theta.
//
// 20140412  Michael Kelsey
// 20140425  Add calculations of E_HV and E_mass (using mass tensor)
// 20140502  Replace NistManager with single-element material creation;
//		avoids weird segfault from NistManager deletion.

#include "G4LatticeLogical.hh"
#include "G4LatticeManager.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <iostream>
using namespace std;

G4bool use_valley = false;	// Set false to use HV wavevectors

int main(/*int argc, char* argv[]*/) {
  G4Material* ge = new G4Material("Ge", 32., 72.630*g/mole, 5.323*g/cm3,
				  kStateSolid);
  G4LatticeLogical* lattice =
    G4LatticeManager::GetLatticeManager()->LoadLattice(ge, "Ge");
  G4double me_HV = lattice->GetElectronMass();

  // Electron wavevector corresponding to phonon speed in lattice
  G4double ksound = lattice->GetSoundSpeed() * me_HV / hbar_Planck;

  // Buffers to construct wavevectors and momenta
  G4ThreeVector k_HV, q_HV, k_recoil;
  G4ThreeVector p, qMom, p_recoil, delta_p;
  G4ThreeVector p_valley, psq;
  G4ThreeVector k_valley, q_v;

  const G4RotationMatrix& mass = lattice->GetMassTensor();

  // Column headings for analysis
  cout << "kmag\tk/ks\tthetaMax\ttheta\tq_HV\tk_rec\tp_in\tp_phon\tp_recoil"
       << "\tEtrack\tErecoil\tE_NIEL\tE_HV\tE_mass" << endl;

  // Only "supersonic" electrons generate Luke phonons
  G4double dk = 0.1*ksound;
  for (G4double kmag=ksound+dk; kmag<=50*ksound; kmag+=dk) {
    k_HV.setRThetaPhi(kmag, 0., 0.);		// Align wavevector on valley
    k_valley = lattice->MapK_HVtoK_valley(1, k_HV);
    p = lattice->MapK_HVtoP(1, k_HV);

    G4double thMax = acos(ksound/kmag);		// Maximum emission angle
    for (G4double th=0; th<=thMax; th+=0.1) {
      G4double q = 2.*(kmag*cos(th) - ksound);  // Phonon wavevector
      q_v.setRThetaPhi(q, th, 0.);

      if (use_valley) {
	qMom = lattice->MapK_valleyToP(1, q_v);
	k_recoil = k_valley - q_v;
	p_recoil = lattice->MapK_valleyToP(1, k_recoil);
      } else {
	qMom = lattice->MapK_HVtoP(1, q_v);
	k_recoil = k_HV - q_v;
	p_recoil = lattice->MapK_HVtoP(1, k_recoil);
      }

      delta_p = p - qMom;

      G4double Etrack = 0.5*p.mag2()/(me_HV*c_squared);
      G4double Erecoil = 0.5*p_recoil.mag2()/(me_HV*c_squared);

      // Kinetic energy in HV space
      G4ThreeVector p_HV = lattice->MapPtoK_HV(1, p_recoil) * hbar_Planck;
      G4double E_HV = 0.5 * p_HV.mag2() / me_HV;

      // Kinetic energy using mass tensor on components
      p_valley = lattice->GetValley(1) * p_recoil;
      psq.set(p_valley.x()*p_valley.x(), p_valley.y()*p_valley.y(),
	      p_valley.z()*p_valley.z());
      G4double E_mass =
	0.5 * (psq.x()/mass.xx() + psq.y()/mass.yy() + psq.z()/mass.zz()) / c_squared;

      cout << kmag << "\t" << kmag/ksound << "\t" << thMax << "\t" << th
	   << "\t" << q << "\t" << k_recoil.mag() << "\t" << p.mag()
	   << "\t" << qMom.mag() << "\t" << delta_p.mag() //***p_recoil.mag()
	   << "\t" << Etrack << "\t" << Erecoil << "\t" << Etrack-Erecoil
	   << "\t" << E_HV << "\t" << E_mass << endl;
    }	// for (G4double th...
  }	// for (G4double kmag...
}

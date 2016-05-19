/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// latticeVecs.cc       Show valley directions for different lattice
//                      orientations.
//
// 20160519  Michael Kelsey

#include "G4Box.hh"
#include "G4LatticeLogical.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <iostream>
#include <iomanip>
using namespace std;


void showLatticeFrame(G4LatticePhysical* lat, G4ThreeVector vec) {
  G4cout << " Vector " << vec << " in lattice frame "
	 << lat->RotateToLattice(vec) << G4endl;
}

int main(int argc, char* argv[]) {
  // Materials available
  G4Material* ge = new G4Material("Ge", 32., 72.630*g/mole, 5.323*g/cm3,
                                  kStateSolid);
  G4Material* si = new G4Material("Si", 14., 28.085*g/mole, 2.329*g/cm3,
				  kStateSolid);

  // Need physical volumes in order to make physical lattice
  G4Material* vac = new G4Material("Vacuum",1.,1*g/mole,1e-20*g/cm3,kStateGas);
  G4VSolid* worldS = new G4Box("World", 5*cm, 5*cm, 5*cm);
  G4LogicalVolume* world = new G4LogicalVolume(worldS, vac, "World");

  G4VSolid* crystal = new G4Box("Crystal", 1*cm, 1*cm, 1*cm);
  G4LogicalVolume* crystalGe = new G4LogicalVolume(crystal, ge, "GeCrystal");
  G4LogicalVolume* crystalSi = new G4LogicalVolume(crystal, si, "SiCrystal");

  G4VPhysicalVolume* pvGe =
    new G4PVPlacement(0, G4ThreeVector(-2*cm,0.,0.), crystalGe,
		      crystalGe->GetName(), world, 0, 0, false);

  G4VPhysicalVolume* pvSi =
    new G4PVPlacement(0, G4ThreeVector(2*cm,0.,0.), crystalSi,
		      crystalSi->GetName(), world, 0, 0, false);

  // Construct physical lattices for each crystal
  G4LatticePhysical* geLat =
    G4LatticeManager::GetLatticeManager()->LoadLattice(pvGe, "Ge");
  G4cout << "Germainum lattice:\n" << *geLat << G4endl;

  /***
  G4LatticePhysical* siLat =
    G4LatticeManager::GetLatticeManager()->LoadLattice(pvSi, "Si");
  siLat->SetVerboseLevel(1);
  G4cout << "Silicon lattice:\n" << *siLat << G4endl;
  ***/

  geLat->SetVerboseLevel(1);

  // Report what should be null rotations
  G4ThreeVector diag(1.,1.,1.);
  showLatticeFrame(geLat, CLHEP::HepXHat);
  showLatticeFrame(geLat, CLHEP::HepYHat);
  showLatticeFrame(geLat, CLHEP::HepZHat);
  showLatticeFrame(geLat, diag);

  // Put lattice into (100) orientation
  geLat->SetMillerOrientation(1,0,0);
  showLatticeFrame(geLat, CLHEP::HepXHat);
  showLatticeFrame(geLat, CLHEP::HepYHat);
  showLatticeFrame(geLat, CLHEP::HepZHat);
  showLatticeFrame(geLat, diag);

  // Put lattice into (001) orientation
  geLat->SetMillerOrientation(0,0,1);
  showLatticeFrame(geLat, CLHEP::HepXHat);
  showLatticeFrame(geLat, CLHEP::HepYHat);
  showLatticeFrame(geLat, CLHEP::HepZHat);
  showLatticeFrame(geLat, diag);

  // Put lattice into (111) orientation
  geLat->SetMillerOrientation(1,1,1);
  showLatticeFrame(geLat, CLHEP::HepXHat);
  showLatticeFrame(geLat, CLHEP::HepYHat);
  showLatticeFrame(geLat, CLHEP::HepZHat);
  showLatticeFrame(geLat, diag);
}

// testHVtransform: Exercise transform into and out of Herring-Vogt frame
//
// Usage: testHVtransform [-v N]
//
// Options: -v N	Set vebosity to N: 1 = print errors, 2 = print all
//
// 20201127  Michael Kelsey

#include "globals.hh"
#include "G4CMPConfigManager.hh"
#include "G4LatticeLogical.hh"
#include "G4LatticeManager.hh"
#include "G4Material.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"


// Flag to print out all calculations

namespace { G4int verbose = 0; }

// Construct or return lattice for transforms

const G4LatticeLogical* getLattice() {
  static G4LatticeLogical* theLattice = 0;
  if (!theLattice) {
    G4Material* mat = new G4Material("Ge", 32., 72.63, 5.323*g/cm3);
    theLattice = G4LatticeManager::Instance()->LoadLattice(mat, "Ge");
  }

  return theLattice;
}

// Exercise HV transforms on unit vector in lattice frame

G4bool bigDiff(G4double diffval) { return fabs(diffval)>1e-6; }
G4bool showDiff(G4double diffval) {
  return (verbose && (verbose>1 || bigDiff(diffval)));
}

G4int testHVstep(G4double costh, G4double phi) {
  const G4LatticeLogical* lat = getLattice();

  G4ThreeVector mom; 
  mom.setRThetaPhi(10.*eV, acos(costh), phi);

  if (verbose>1) {
    G4cout << "\ncos(theta) " << costh << " phi " 
	   << phi/deg << " deg " << mom/eV << " eV" << G4endl;
  }

  G4int ndiff = 0;

  // Momentum/velocity relationship
  G4ThreeVector vel = lat->MapPtoV_el(1, mom);
  G4ThreeVector vP = lat->MapV_elToP(1, vel);
  G4double pvDiff = (mom-vP).mag();
  if (bigDiff(pvDiff)) ndiff++;

  if (showDiff(pvDiff)) {
    G4cout << " MapPtoV_el " << vel/(km/s) << " km/s and back " << vP/eV
	   << " (" << pvDiff/eV << ")" << G4endl;
  }

  // Momentum/velocity to HV wavevector
  G4ThreeVector pKHV = lat->MapPtoK_HV(1, mom);
  G4ThreeVector kHVp = lat->MapK_HVtoP(1, pKHV);
  G4double pkHVDiff = (mom-kHVp).mag();
  if (bigDiff(pkHVDiff)) ndiff++;

  if (showDiff(pkHVDiff)) {
    G4cout << " MapPtoK_HV " << pKHV*um << "/um and back " << kHVp/eV
	   << " (" << pkHVDiff/eV << ")" << G4endl;
  }

  G4ThreeVector vKHV = lat->MapV_elToK_HV(1, vel);
  G4double kvkDiff = (pKHV-vKHV).mag();

  if (showDiff(kvkDiff)) {
    G4cout << " MapV_elToK_HV " << vKHV*um << "/um (" << kvkDiff*um << ")"
	   << G4endl;
  }

  // Wavevectors
  G4ThreeVector pK = lat->MapPtoK_valley(1, mom);
  G4ThreeVector kP = lat->MapK_valleyToP(1, pK);
  G4double pkDiff = (mom-kP).mag();
  if (bigDiff(pkDiff)) ndiff++;

  if (showDiff(pkDiff)) {
    G4cout << " MapPtoK_valley " << pK*um << "/um and back " << kP/eV
	   << " (" << pkDiff/eV << ")" << G4endl;
  }

  G4ThreeVector kHV_K = lat->MapK_HVtoK_valley(1, pKHV);
  G4double kkDiff = (pK-kHV_K).mag();
  if (bigDiff(kkDiff)) ndiff++;

  if (showDiff(kkDiff)) {
    G4cout << " MapK_HVtoK_valley " << kHV_K*um << "/um (" << kkDiff*m << ")"
	   << G4endl;
  }

  // Momentum or velocity to kinetic energy
  G4double pE = lat->MapPtoEkin(1, mom);
  G4double vE = lat->MapV_elToEkin(1, vel);
  if (bigDiff(pE-vE)) ndiff++;

  if (showDiff(pE-vE)) {
    G4cout << " MapPtoEkin    " << pE/eV << " eV MapV_eltoEkin " << vE/eV
	   << " eV (" << (pE-vE)/eV << ")" << G4endl;
  }

  return ndiff;
}

// Scan polar angle in rotation phi around Z axis

G4int testHVpolar(G4double costh) {
  G4int ndiff = 0;
  for (G4double phi=0.; phi<360.; phi+=30.) {
    G4int ntry = testHVstep(costh, phi*deg);
    if (ntry>0 && verbose) 
      G4cerr << "Got " << ntry << " discrepancies" << G4endl;

    ndiff += ntry;
  }

  return ndiff;
}


// MAIN PROGRAM

int main(int argc, char* argv[]) {
  // Get optional command line argument for verbosity
  if (argc > 2 && strcmp(argv[1],"-v")==0) {
    verbose = atoi(argv[2]);
  }

  G4CMPConfigManager::SetVerboseLevel(verbose);

  if (verbose) {
    G4cout << "testHVtransform:  Exercise Herring-Vogt transforms.  Uses\n"
	   << "vectors in lattice frame to transform into, then back out\n"
	   << "of, Herring-Vogt frame to confirm invertibility." << G4endl;
  }

  // Take uniform steps in cos(theta)
  G4int ndiff = 0;
  for (G4double costh=1.; costh>=-1; costh-=0.2) {
    ndiff += testHVpolar(costh);
  }

  return ndiff;		// 0 is successful validation
}

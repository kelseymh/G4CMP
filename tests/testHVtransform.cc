// testHVtransform: Exercise transform into and out of Herring-Vogt frame
//
// Usage: testHVtransform [-v N]
//
// Options: -v N	Set vebosity to N: 1 = print errors, 2 = print all
//
// 20201127  Michael Kelsey
// 20210921  Look at step-wise transforms (global-local-lattice-valley-HV)
//		as used, for instance, in G4CMPEqEMField.

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


// Scan cos(theta) angle in steps relative to valley axis

G4int testHVangle(G4double costh) {
  const G4LatticeLogical* theLattice = getLattice();

  G4int ndiff=0;

  G4ThreeVector valmom(costh, 0., 1-costh*costh);	// X-Z plane
  G4ThreeVector valHV = theLattice->GetSqrtInvTensor()*valmom;

  if (verbose) {
    G4ThreeVector xhat(1,0,0);
    G4cout << " cos(theta) " << costh << " -> cos(thHV) "
	   << valHV.unit().dot(xhat) << " mag " << valHV << G4endl;
  }

  G4double magdiff = valHV.mag() - valmom.mag();
  if (bigDiff(magdiff)) ndiff++;
  if (showDiff(magdiff)) {
    G4cout << " unit vector at " << costh << " has HV magnitude "
	   << valHV.mag() << G4endl;
  }

  return ndiff;
}


// Apply lattice-to-valley transforms individually, compare with expectations

G4int testValleyAxis(G4int iv) {
  const G4LatticeLogical* lat = getLattice();	// Germanium
  if (iv < 0 || iv >= (G4int)lat->NumberOfValleys()) return 1000;

  G4int ndiff=0;

  G4ThreeVector xhat(1,0,0);	    // Valley axis specifies local X direction

  G4ThreeVector axis = lat->GetValleyAxis(iv);
  if (verbose) G4cout << "Valley " << iv << " axis " << axis << G4endl;
  
  G4ThreeVector vhat = lat->GetValley(iv)*axis;	// SHOULD EQUAL xhat
  G4double xvdot = xhat.dot(vhat);
  if (bigDiff(xvdot)) ndiff++;
  if (showDiff(xvdot)) {
    G4cout << "GetValley(iv)*axis " << vhat << " vs xhat = " << xvdot << G4endl;
  }
  
  G4ThreeVector vloc = lat->GetValleyInv(iv)*xhat;	// SHOULD EQUAL axis
  G4double avdot = axis.dot(vloc);
  if (bigDiff(avdot)) ndiff++;
  if (showDiff(avdot)) {
    G4cout << "GetValleyInv(iv)*xhat " << vloc << " vs axis = " << avdot << G4endl;
  }

  return ndiff;
}

// Test off-axis transforms out of valley frame (see G4CMP-279)

G4int testValleyOffAxis(G4int iv) {
  const G4LatticeLogical* lat = getLattice();	// Germanium
  if (iv < 0 || iv >= (G4int)lat->NumberOfValleys()) return 1000;

  G4int ndiff = 0;

  G4ThreeVector axis = lat->GetValleyAxis(iv);

  // Start in valley frame, along Z with some axial component
  G4ThreeVector xzvalley(0.2,0.,0.98); xzvalley.setMag(1.);
  G4double xzvAngle = atan2(xzvalley.z(),xzvalley.x());
  if (verbose) {
    G4cout << "Valley " << iv << " xzv " << xzvalley << " angle "
	   << xzvAngle/deg << " deg to axis" << G4endl;
  }

  G4ThreeVector xzvLat = lat->GetValleyInv(iv)*xzvalley;
  G4double xzlAngle = acos(axis.dot(xzvLat));
  if (verbose) {
    G4cout << " in lattice becomes " << xzvLat << " angle " << xzlAngle/deg
	   << " deg to axis " << G4endl;
  }

  if (bigDiff(xzlAngle-xzvAngle)) ndiff++;

  return ndiff;
}

// Apply valley transform to Zhat (e.g. E-field) show result

G4int testZhatValley(G4int iv) {
  const G4LatticeLogical* lat = getLattice();	// Germanium
  if (iv < 0 || iv >= (G4int)lat->NumberOfValleys()) return 1000;

  G4int ndiff = 0;

  G4ThreeVector zhat(0,0,1);		// Typical E-field in experiments
  G4ThreeVector zvalley = lat->GetValley(iv)*zhat;

  if (verbose)
    G4cout << "Valley " << iv << " zhat looks like " << zvalley << G4endl;

  // Inverse transform should recover original zhat
  G4ThreeVector zloc = lat->GetValleyInv(iv)*zvalley;
  G4double zdiff = zhat.dot(zloc);
  ndiff += bigDiff(zdiff);

  if (showDiff(zdiff)) {
    G4cout << " inverse transform gives " << zloc << " vs zhat " << zdiff
	   << G4endl;
  }

  return ndiff;
}

// Apply valley transform to Zhat along with 1/minv momentum transform

G4int testZminvValley(G4int iv) {
  const G4LatticeLogical* lat = getLattice();	// Germanium
  if (iv < 0 || iv >= (G4int)lat->NumberOfValleys()) return 1000;

  G4int ndiff = 0;

  G4ThreeVector zhat(0,0,1);		// Typical E-field in experiments
  G4ThreeVector zvalley = lat->GetValley(iv)*zhat;
  G4ThreeVector vvalley = lat->GetMInvTensor()*zvalley * lat->GetElectronMass();
  G4ThreeVector vHV = lat->GetSqrtInvTensor()*zvalley;

  if (verbose) {
    G4cout << "Valley " << iv << " zhat/M     " << vvalley << G4endl
	   << "Valley " << iv << " zhat/sqrtM " << vHV << G4endl;
  }

  // Rotate transformed vector back into lattice frame
  G4ThreeVector vlat = lat->GetValleyInv(iv)*vHV;
  if (verbose)
    G4cout << " rotated vHV back to lattice, vlat " << vlat << G4endl;

  G4double ediff = vlat.mag() - zhat.mag();
  if (zhat.dot(vlat) < 0. || bigDiff(ediff)) ndiff++;
  if (showDiff(ediff))
    G4cout << " Vector magnitude changed by " << ediff << G4endl;

  return ndiff;
}

// Apply valley transform to Yhat (e.g. E-field in tilted frame) show result

G4int testYhatValley(G4int iv) {
  const G4LatticeLogical* lat = getLattice();	// Germanium
  if (iv < 0 || iv >= (G4int)lat->NumberOfValleys()) return 1000;

  G4int ndiff = 0;

  G4ThreeVector yhat(0,1,0);		// Typical E-field in experiments
  G4ThreeVector yvalley = lat->GetValley(iv)*yhat;

  if (verbose)
    G4cout << "Valley " << iv << " yhat looks like " << yvalley << G4endl;

  // Inverse transform should recover original yhat
  G4ThreeVector yloc = lat->GetValleyInv(iv)*yvalley;
  G4double ydiff = yhat.dot(yloc);
  ndiff += bigDiff(ydiff);

  if (showDiff(ydiff)) {
    G4cout << " inverse transform gives " << yloc << " vs yhat " << ydiff
	   << G4endl;
  }

  return ndiff;
}

// Apply valley transform to Yhat along with 1/minv momentum transform

G4int testYminvValley(G4int iv) {
  const G4LatticeLogical* lat = getLattice();	// Germanium
  if (iv < 0 || iv >= (G4int)lat->NumberOfValleys()) return 1000;

  G4int ndiff = 0;

  G4ThreeVector yhat(0,1,0);		// Typical E-field in experiments
  G4ThreeVector yvalley = lat->GetValley(iv)*yhat;
  G4ThreeVector vvalley = lat->GetMInvTensor()*yvalley * lat->GetElectronMass();
  G4ThreeVector vHV = lat->GetSqrtInvTensor()*yvalley;

  if (verbose) {
    G4cout << "Valley " << iv << " yhat/M     " << vvalley << G4endl
	   << "Valley " << iv << " yhat/sqrtM " << vHV << G4endl;
  }

  // Rotate transformed vector back into lattice frame
  G4ThreeVector vlat = lat->GetValleyInv(iv)*vHV;
  if (verbose)
    G4cout << " rotated vHV back to lattice, vlat " << vlat << G4endl;

  G4double ediff = vlat.mag() - yhat.mag();
  if (yhat.dot(vlat) < 0. || bigDiff(ediff)) ndiff++;
  if (showDiff(ediff))
    G4cout << " Vector magnitude changed by " << ediff << G4endl;

  return ndiff;
}


// Test various transforms involving individual valleys

G4int testValleys() {
  const G4LatticeLogical* lat = getLattice();	// Germanium

  // Germanium valley axes should all be along (+-1, +-1, +-1)
  G4int ndiff=0;
  for (size_t iv=0; iv<lat->NumberOfValleys(); iv++) {
    if (verbose) G4cout << G4endl;		// Spacing between tests
    ndiff += testValleyAxis(iv);
    ndiff += testValleyOffAxis(iv);
    ndiff += testZhatValley(iv);
    ndiff += testZminvValley(iv);
    ndiff += testYhatValley(iv);
    ndiff += testYminvValley(iv);
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
  for (G4double costh=1.; costh>=-1; costh-=0.1) {
    ndiff += testHVangle(costh);
  }

  // Check valley axes and to-from valley transforms
  ndiff += testValleys();

  return ndiff;		// 0 is successful validation
}

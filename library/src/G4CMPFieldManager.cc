// $Id$
//
// 20140328  Save field locally, get local/global transform from track to
//	     pass to field and to chord finder for electrons.
// 20140329  Pass G4CMP field pointer, which handles local/global transform
// 20140331  Use G4CMPEqEMField for everything, now handles holes; don't
//	     need lattice locally; get physical lattice track by track

#include "G4CMPFieldManager.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPEqEMField.hh"
#include "G4CMPLocalElectroMagField.hh"
#include "G4CMPValleyTrackMap.hh"
#include "G4ChordFinder.hh"
#include "G4ClassicalRK4.hh"
#include "G4ElectroMagneticField.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"


// Constructors and destructor

G4CMPFieldManager::G4CMPFieldManager(G4ElectroMagneticField *detectorField)
  : G4FieldManager(detectorField),
    myDetectorField(new G4CMPLocalElectroMagField(detectorField)),
    stepperVars(8), stepperLength(1e-9*mm) {
  CreateTransport();
}

G4CMPFieldManager::G4CMPFieldManager(G4CMPLocalElectroMagField *detectorField)
  : G4FieldManager(detectorField), myDetectorField(detectorField),
    stepperVars(8), stepperLength(1e-9*mm) {
  CreateTransport();
}

void G4CMPFieldManager::CreateTransport() {
  theEqMotion    = new G4CMPEqEMField(myDetectorField);
  theStepper     = new G4ClassicalRK4(theEqMotion, stepperVars);
  theDriver      = new G4MagInt_Driver(stepperLength, theStepper, stepperVars);
  theChordFinder = new G4ChordFinder(theDriver);
  SetChordFinder(theChordFinder);
}

G4CMPFieldManager::~G4CMPFieldManager() {
  delete myDetectorField;   myDetectorField=0;
  delete theEqMotion;       theEqMotion=0;
  delete theStepper;        theStepper=0;
  delete theDriver;         theDriver=0;
  delete theChordFinder;    theChordFinder=0;
}


void G4CMPFieldManager::ConfigureForTrack(const G4Track* aTrack) {
  // Configure equation of motion with physical lattice
  const G4LatticePhysical* lat =
    G4LatticeManager::GetLatticeManager()->GetLattice(aTrack->GetVolume());
  G4bool newLat = theEqMotion->ChangeLattice(lat);

  // Extract local/global transform from track only if volume has changed
  // This avoids creating and copying transforms on every step
  if (newLat) {
    const G4RotationMatrix* rot = aTrack->GetTouchable()->GetRotation();
    const G4ThreeVector& trans  = aTrack->GetTouchable()->GetTranslation();
    G4AffineTransform localToGlobal(rot, trans);

    myDetectorField->SetTransforms(localToGlobal);
    theEqMotion->SetTransforms(localToGlobal);
  }

  // Configure electric field with valleys for either electrons or holes
  if (aTrack->GetDefinition() == G4CMPDriftElectron::Definition()) {
    G4int iv = G4CMPValleyTrackMap::GetInstance()->GetValley(aTrack);
    theEqMotion->SetValley(iv);
  } else {
    theEqMotion->SetNoValley();
  }
}



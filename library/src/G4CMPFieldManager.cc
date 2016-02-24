/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
//
// 20140328  Save field locally, get local/global transform from track to
//	     pass to field and to chord finder for electrons.
// 20140329  Pass G4CMP field pointer, which handles local/global transform
// 20140331  Use G4CMPEqEMField for everything, now handles holes; don't
//	     need lattice locally; get physical lattice track by track
// 20150122  Use verboseLevel instead of compiler flag for debugging
// 20150528  Pass verbosity through to field computation classes

#include "G4CMPFieldManager.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPEqEMField.hh"
#include "G4CMPLocalElectroMagField.hh"
#include "G4CMPTrackInformation.hh"
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
  : G4FieldManager(new G4CMPLocalElectroMagField(detectorField)),
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

  theDriver->SetVerboseLevel(G4CMPConfigManager::GetVerboseLevel());
  theChordFinder->SetVerbose(G4CMPConfigManager::GetVerboseLevel());
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

    if (G4CMPConfigManager::GetVerboseLevel() > 1) {
      G4cout << "G4CMPFieldManager::ConfigureForTrack with translation "
	     << trans << " rotation " << *rot << G4endl;
    }

    myDetectorField->SetTransforms(localToGlobal);
    theEqMotion->SetTransforms(localToGlobal);
  }

  // Configure electric field with valleys for either electrons or holes
  if (aTrack->GetDefinition() == G4CMPDriftElectron::Definition()) {
    G4int modelID = G4PhysicsModelCatalog::GetIndex("G4CMP process");
    if (modelID < 0) {
      G4Exception("G4CMPFieldManager::ConfigureForTrack","Electron001",
      EventMustBeAborted, "Track is electron, but has no G4CMP Aux. Info.");
    }
    G4int iv =
      static_cast<G4CMPTrackInformation*>(
        aTrack->GetAuxiliaryTrackInformation(modelID)
                                         )->GetValleyIndex();
    SetElectronValleyForTrack(iv);
  } else {
    SetElectronValleyForTrack(-1);
  }
}

void G4CMPFieldManager::SetElectronValleyForTrack(G4int valley) {
  if(valley >= 0 && valley <= 3) {
    theEqMotion->SetValley(valley);
  } else {
    theEqMotion->SetNoValley();
  }
}



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
// 20161114  Use new G4CMPDriftTrackInfo
// 20170801  Count consecutive null lattice pointers for reflection steps
// 20180201  Add G4MagIntegratorDriver.hh, needed with Geant4 10.4.
// 20180319  Don't delete theDriver; done by G4ChordFinder.
// 20200213  In ConfigureForTrack, check if registered field is wrapped in
//		G4CMPLocalEMField; apply wrapping if needed.
// 20200804  Attach local geometry shape to field
// 20210901  Add local verbosity flag for reporting diagnostics, pass through
//		to G4CMPLocalEMField.

#include "G4CMPFieldManager.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4CMPEqEMField.hh"
#include "G4CMPLocalElectroMagField.hh"
#include "G4CMPDriftTrackInfo.hh"
#include "G4CMPTrackUtils.hh"
#include "G4ChordFinder.hh"
#include "G4ClassicalRK4.hh"
#include "G4ElectroMagneticField.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalVolume.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"
#include "G4VTouchable.hh"


// Constructors and destructor

G4CMPFieldManager::G4CMPFieldManager(G4ElectroMagneticField *detectorField,
				     G4int vb)
  : G4FieldManager(new G4CMPLocalElectroMagField(detectorField)),
    verboseLevel(vb==0?G4CMPConfigManager::GetVerboseLevel():vb),
    myDetectorField(0), stepperVars(8), stepperLength(1e-9*mm),
    latticeNulls(0), maxLatticeNulls(3) {
  if (verboseLevel)
    G4cout << "G4CMPFieldManager wrapped global field in LocalEMField." << G4endl;

  // Same pointer, but non-const for use in ConfigureForTrack()
  G4Field* baseField = const_cast<G4Field*>(GetDetectorField());
  myDetectorField = dynamic_cast<G4CMPLocalElectroMagField*>(baseField);
  if (myDetectorField) myDetectorField->SetVerboseLevel(verboseLevel);

  CreateTransport();
}

G4CMPFieldManager::G4CMPFieldManager(G4CMPLocalElectroMagField *detectorField,
				     G4int vb)
  : G4FieldManager(detectorField),
    verboseLevel(vb==0?G4CMPConfigManager::GetVerboseLevel():vb),
    myDetectorField(detectorField), stepperVars(8), stepperLength(1e-9*mm),
    latticeNulls(0), maxLatticeNulls(3) {
  if (verboseLevel)
    G4cout << "G4CMPFieldManager provided with wrapped LocalEMField." << G4endl;

  if (myDetectorField) myDetectorField->SetVerboseLevel(verboseLevel);
  CreateTransport();
}

G4CMPFieldManager::~G4CMPFieldManager() {
  delete theEqMotion;       theEqMotion=0;
  delete theStepper;        theStepper=0;
  delete theChordFinder;    theChordFinder=0;
}


// Create all of the pieces for transport through electric field

void G4CMPFieldManager::CreateTransport() {
  theEqMotion    = new G4CMPEqEMField(myDetectorField);
  theStepper     = new G4ClassicalRK4(theEqMotion, stepperVars);
  theDriver      = new G4MagInt_Driver(stepperLength, theStepper, stepperVars);
  theChordFinder = new G4ChordFinder(theDriver);
  SetChordFinder(theChordFinder);

  theEqMotion->SetVerboseLevel(verboseLevel);
  theDriver->SetVerboseLevel(verboseLevel);
  theChordFinder->SetVerbose(verboseLevel);
}


// Run-time Configuration

void G4CMPFieldManager::ConfigureForTrack(const G4Track* aTrack) {
  if (verboseLevel) {
    G4cout << "G4CMPFieldManager::ConfigureForTrack "
	   << aTrack->GetTrackID() << "/" << aTrack->GetCurrentStepNumber()
	   << " @ " << aTrack->GetPosition() << " in "
	   << aTrack->GetVolume()->GetName() << G4endl;
  }

  // Ensure that field is properly wrapped for global/local coordinates
  if (!dynamic_cast<const G4CMPLocalElectroMagField*>(GetDetectorField())) {
    if (verboseLevel) {
      G4cout << " Registered field not local.  Wrapping in G4CMPLocalEMField."
	     << G4endl;
    }

    const G4ElectroMagneticField* baseField =
      dynamic_cast<const G4ElectroMagneticField*>(GetDetectorField());
    if (!baseField) {
      G4ExceptionDescription msg;
      msg << "Field attached to volume " << aTrack->GetVolume()->GetName()
	  << " not G4ElectroMagneticField.";

      G4Exception("G4CMPFieldManager::ConfigureForTrack", "FieldMan003",
		  FatalException, msg);
      return;
    }

    myDetectorField = new G4CMPLocalElectroMagField(baseField);
    myDetectorField->SetVerboseLevel(verboseLevel);

    ChangeDetectorField(myDetectorField);
  }

  // Configure equation of motion with physical lattice
  const G4LatticePhysical* lat =
    G4LatticeManager::GetLatticeManager()->GetLattice(aTrack->GetVolume());

  // If track is outside valid volume, count attemps to look for reflections
  if (lat) latticeNulls = 0;
  else {
    latticeNulls++;
    if (latticeNulls > maxLatticeNulls) {
      G4ExceptionDescription msg;
      msg << "No lattice available for volume "
	  << aTrack->GetVolume()->GetName() << " after " << latticeNulls
	  << " steps.";

      G4Exception("G4CMPFieldManager::ConfigureForTrack", "FieldMan002",
		  EventMustBeAborted, msg);
    }

    theEqMotion->SetNoValley();
    return;
  }

  // Hack around boundary issues; don't store or change vol if null lattice!
  G4bool newLat = lat ? theEqMotion->ChangeLattice(lat) : false;

  // Extract local/global transform from track only if volume has changed
  // This avoids creating and copying transforms on every step
  if (newLat) {
    const G4RotationMatrix* rot = aTrack->GetTouchable()->GetRotation();
    const G4ThreeVector& trans  = aTrack->GetTouchable()->GetTranslation();
    G4AffineTransform localToGlobal(rot, trans);

    if (verboseLevel > 1) {
      G4cout << " volume " << aTrack->GetVolume()->GetName() << " has"
	     << " trans " << trans << " rot " << rot->delta()/deg
	     << " deg about " << rot->axis() << G4endl;
    }

    myDetectorField->SetGeometry(aTrack->GetTouchable()->GetVolume()->
				 GetLogicalVolume()->GetSolid());
    myDetectorField->SetTransforms(localToGlobal);
    theEqMotion->SetTransforms(localToGlobal);
  }

  G4int iv = -1;
  if (lat && (aTrack->GetDefinition() == G4CMPDriftElectron::Definition() ||
	      aTrack->GetDefinition() == G4CMPDriftHole::Definition())) {
    iv = G4CMP::GetTrackInfo<G4CMPDriftTrackInfo>(*aTrack)->ValleyIndex();
    SetChargeValleyForTrack(lat, iv);
  } else {
    theEqMotion->SetNoValley();
  }
}

void G4CMPFieldManager::SetChargeValleyForTrack(const G4LatticePhysical* lat,
                                                G4int valley) {
  if (valley < -1 || valley > static_cast<G4int>(lat->NumberOfValleys() - 1)) {
    G4Exception("G4CMPFieldManager::SetChargeValleyForTrack", "FieldMan001",
                EventMustBeAborted,
                "Valley index is not valid for the current lattice.");
    return;
  }

  if(valley >= 0) {
    theEqMotion->SetValley(valley);
  } else {
    theEqMotion->SetNoValley();
  }
}



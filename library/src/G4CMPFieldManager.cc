#include "G4CMPFieldManager.hh"
#include "G4CMPEqEMField.hh"
#include "G4AffineTransform.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPValleyTrackMap.hh"
#include "G4ChordFinder.hh"
#include "G4ClassicalRK4.hh"
#include "G4ElectricField.hh"
#include "G4EqMagElectricField.hh"
#include "G4LatticeLogical.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"


// Constructor and destructor

G4CMPFieldManager::G4CMPFieldManager(G4ElectricField *detectorField,
				     const G4LatticeLogical* lattice)
  : G4FieldManager(detectorField), stepperVars(8), stepperLength(1e-9*mm),
    EqNormal(0), normalStepper(0), normalDriver(0), normalChordFinder(0),
    theLattice(lattice), nValleys(lattice?lattice->NumberOfValleys():0) {
  // Set up default field for holes, non-valley charged particles
  EqNormal = new G4EqMagElectricField(detectorField);
  normalStepper  = new G4ClassicalRK4(EqNormal, stepperVars);
  normalDriver  = new G4MagInt_Driver(1e-9*mm, normalStepper, stepperVars);
  normalChordFinder  = new G4ChordFinder(normalDriver);

  // Set up fielding handling for valleys, if lattice was provided
  for (size_t iv=0; iv<nValleys; iv++) {
    EqValley.push_back(new G4CMPEqEMField(detectorField, theLattice->GetValley(iv)));
    valleyStepper.push_back(new G4ClassicalRK4(EqValley.back(), stepperVars));
    valleyDriver.push_back(new G4MagInt_Driver(stepperLength, valleyStepper.back(), stepperVars));
    valleyChordFinder.push_back(new G4ChordFinder(valleyDriver.back()));
  }
}

G4CMPFieldManager::~G4CMPFieldManager() {
  delete EqNormal;          EqNormal=0;
  delete normalStepper;     normalStepper=0;
  delete normalDriver;      normalDriver=0;
  delete normalChordFinder; normalChordFinder=0;

  for (size_t iv=0; iv<nValleys; iv++) {
    delete EqValley[iv];          EqValley[iv]=0;
    delete valleyStepper[iv];     valleyStepper[iv]=0;
    delete valleyDriver[iv];      valleyDriver[iv]=0;
    delete valleyChordFinder[iv]; valleyChordFinder[iv]=0;
  }
}


void G4CMPFieldManager::ConfigureForTrack(const G4Track* aTrack) {
  if (aTrack->GetDefinition() == G4CMPDriftElectron::Definition()) {
    G4int iv = G4CMPValleyTrackMap::GetInstance()->GetValley(aTrack);
    if (iv>=0 && iv<(G4int)nValleys) {
      SetChordFinder(valleyChordFinder[iv]);
    } else {
      G4cerr << "ERROR: G4CMPFieldManager: invalid valley " << iv << G4endl;
      SetChordFinder(normalChordFinder);
    }
  } else {
    SetChordFinder(normalChordFinder);
  }
}



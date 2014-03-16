#include "G4CMPFieldManager.hh"
#include "G4CMPEqEMField.hh"
#include "G4AffineTransform.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPValleyTrackMap.hh"
#include "G4ChordFinder.hh"
#include "G4ClassicalRK4.hh"
#include "G4ElectricField.hh"
#include "G4EqMagElectricField.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"


// Constructor and destructor

G4CMPFieldManager::G4CMPFieldManager(G4ElectricField *detectorField)
  : G4FieldManager(detectorField) {
  G4RotationMatrix Mvalley;

  EqNormal = new G4EqMagElectricField(detectorField);

  Mvalley.set(-pi/4,-pi/4, pi/4);
  EqValley1 = new G4CMPEqEMField(detectorField, G4AffineTransform(Mvalley));

  Mvalley.set( pi/4,-pi/4,-pi/4);
  EqValley2 = new G4CMPEqEMField(detectorField, G4AffineTransform(Mvalley));

  Mvalley.set(-pi/4, pi/4, pi/4);
  EqValley3 = new G4CMPEqEMField(detectorField, G4AffineTransform(Mvalley));

  Mvalley.set( pi/4, pi/4,-pi/4);
  EqValley4 = new G4CMPEqEMField(detectorField, G4AffineTransform(Mvalley));
  
  const G4int stepperVars = 8;

  normalStepper  = new G4ClassicalRK4(EqNormal, stepperVars);
  valley1Stepper = new G4ClassicalRK4(EqValley1, stepperVars);
  valley2Stepper = new G4ClassicalRK4(EqValley2, stepperVars);
  valley3Stepper = new G4ClassicalRK4(EqValley3, stepperVars);
  valley4Stepper = new G4ClassicalRK4(EqValley4, stepperVars);
  
  normalDriver  = new G4MagInt_Driver(1e-9*mm, normalStepper, stepperVars);
  valley1Driver = new G4MagInt_Driver(1e-9*mm, valley1Stepper, stepperVars);
  valley2Driver = new G4MagInt_Driver(1e-9*mm, valley2Stepper, stepperVars);
  valley3Driver = new G4MagInt_Driver(1e-9*mm, valley3Stepper, stepperVars);
  valley4Driver = new G4MagInt_Driver(1e-9*mm, valley4Stepper, stepperVars);
  
  normalChordFinder  = new G4ChordFinder(normalDriver);
  valley1ChordFinder = new G4ChordFinder(valley1Driver);
  valley2ChordFinder = new G4ChordFinder(valley2Driver);
  valley3ChordFinder = new G4ChordFinder(valley3Driver);
  valley4ChordFinder = new G4ChordFinder(valley4Driver);
}

G4CMPFieldManager::~G4CMPFieldManager() {
  delete EqNormal;
  delete EqValley1;
  delete EqValley2;
  delete EqValley3;
  delete EqValley4;

  delete normalStepper;
  delete valley1Stepper;
  delete valley2Stepper;
  delete valley3Stepper;
  delete valley4Stepper;

  delete normalDriver;
  delete valley1Driver;
  delete valley2Driver;
  delete valley3Driver;
  delete valley4Driver;

  delete normalChordFinder;
  delete valley1ChordFinder;
  delete valley2ChordFinder;
  delete valley3ChordFinder;
  delete valley4ChordFinder;
}


void G4CMPFieldManager::ConfigureForTrack(const G4Track* aTrack){
  if (aTrack->GetDefinition() == G4CMPDriftElectron::Definition()) {
    int valley = G4CMPValleyTrackMap::GetInstance()->GetValley(aTrack);
    switch(valley) {
    case 1: SetChordFinder(valley1ChordFinder); break;
    case 2: SetChordFinder(valley2ChordFinder); break;
    case 3: SetChordFinder(valley3ChordFinder); break;
    case 4: SetChordFinder(valley4ChordFinder); break;
    }
  } else {
    SetChordFinder(normalChordFinder);
  }
}



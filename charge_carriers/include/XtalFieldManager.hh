
#ifndef XtalFieldManager_h
#define XtalFieldManager_h 1

#define PI 3.1415926

#include "globals.hh"
#include "G4FieldManager.hh"
#include "G4ElectricField.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ChordFinder.hh"
#include "G4ClassicalRK4.hh"
#include "G4AffineTransform.hh"
#include "G4RotationMatrix.hh"

#include "EqEMFieldXtal.hh"
#include "G4EqMagElectricField.hh"

class G4Field;
class G4ChordFinder;


class XtalFieldManager : public G4FieldManager {

public:
  G4EqMagElectricField *EqNormal;
  EqEMFieldXtal *EqValley1;
  EqEMFieldXtal *EqValley2;
  EqEMFieldXtal *EqValley3;
  EqEMFieldXtal *EqValley4;

  
  G4MagIntegratorStepper *normalStepper;
  G4MagIntegratorStepper *valley1Stepper;
  G4MagIntegratorStepper *valley2Stepper;
  G4MagIntegratorStepper *valley3Stepper;
  G4MagIntegratorStepper *valley4Stepper;


  G4MagInt_Driver *normalDriver;
  G4MagInt_Driver *valley1Driver;
  G4MagInt_Driver *valley2Driver;
  G4MagInt_Driver *valley3Driver;
  G4MagInt_Driver *valley4Driver;

  G4ChordFinder *normalChordFinder;
  G4ChordFinder *valley1ChordFinder;
  G4ChordFinder *valley2ChordFinder;
  G4ChordFinder *valley3ChordFinder;
  G4ChordFinder *valley4ChordFinder;

public:
  XtalFieldManager(G4ElectricField *detectorField, G4ChordFinder *pChordFinder, G4bool b) : G4FieldManager(detectorField, pChordFinder, b){
    
    EqNormal = new G4EqMagElectricField(detectorField);
    EqValley1 = new EqEMFieldXtal(detectorField);
    EqValley2 = new EqEMFieldXtal(detectorField);
    EqValley3 = new EqEMFieldXtal(detectorField);
    EqValley4 = new EqEMFieldXtal(detectorField);

    EqValley1->SetValleyTransform(G4AffineTransform(G4RotationMatrix(PI/4, 
-PI/4, PI/4)));
    EqValley2->SetValleyTransform(G4AffineTransform(G4RotationMatrix(PI/4, 
-PI/4, -PI/4)));
    EqValley3->SetValleyTransform(G4AffineTransform(G4RotationMatrix(-PI/4, 
PI/4, PI/4)));
    EqValley4->SetValleyTransform(G4AffineTransform(G4RotationMatrix(PI/4, 
PI/4, -PI/4)));
    
    normalStepper = new G4ClassicalRK4(EqNormal, 8);
    valley1Stepper = new G4ClassicalRK4(EqValley1, 8);
    valley2Stepper = new G4ClassicalRK4(EqValley2, 8);
    valley3Stepper = new G4ClassicalRK4(EqValley3, 8);
    valley4Stepper = new G4ClassicalRK4(EqValley4, 8);

    normalDriver = new G4MagInt_Driver(1e-6*mm, 
                                       normalStepper, 
                                       normalStepper->GetNumberOfVariables() );

    valley1Driver = new G4MagInt_Driver(1e-6*mm, 
                                        valley1Stepper, 
                                        valley1Stepper->GetNumberOfVariables() );

    valley2Driver = new G4MagInt_Driver(1e-6*mm, 
                                        valley2Stepper, 
                                        valley2Stepper->GetNumberOfVariables() );

    valley3Driver = new G4MagInt_Driver(1e-6*mm, 
                                        valley3Stepper, 
                                        valley3Stepper->GetNumberOfVariables() );

    valley4Driver = new G4MagInt_Driver(1e-6*mm, 
                                        valley4Stepper, 
                                        valley4Stepper->GetNumberOfVariables() );



    normalChordFinder = new G4ChordFinder(normalDriver);
    valley1ChordFinder = new G4ChordFinder(valley1Driver);
    valley2ChordFinder = new G4ChordFinder(valley2Driver);
    valley3ChordFinder = new G4ChordFinder(valley3Driver);
    valley4ChordFinder = new G4ChordFinder(valley4Driver);

}  
  void ConfigureForTrack(const G4Track* aTrack);

};


#endif

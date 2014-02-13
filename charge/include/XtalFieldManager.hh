
#ifndef XtalFieldManager_h
#define XtalFieldManager_h 1

#include "globals.hh"
#include "G4FieldManager.hh"

class DMCClassicalRK4;
class EqEMFieldXtal;
class G4ChordFinder;
class G4ClassicalRK4;
class G4ElectricField;
class G4EqMagElectricField;
class G4MagInt_Driver;
class G4MagIntegratorStepper;

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
  XtalFieldManager(G4ElectricField *detectorField,
		   G4ChordFinder *pChordFinder, G4bool b);
  ~XtalFieldManager();

  void ConfigureForTrack(const G4Track* aTrack);
};


#endif

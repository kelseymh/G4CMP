// $Id$

#ifndef G4CMPFieldManager_h
#define G4CMPFieldManager_h 1

#include "globals.hh"
#include "G4FieldManager.hh"

class G4CMPEqEMField;
class G4ChordFinder;
class G4ElectricField;
class G4EqMagElectricField;
class G4MagInt_Driver;
class G4MagIntegratorStepper;


class G4CMPFieldManager : public G4FieldManager {
public:
  G4CMPFieldManager(G4ElectricField *detectorField);
  ~G4CMPFieldManager();

  void ConfigureForTrack(const G4Track* aTrack);

private:
  G4EqMagElectricField *EqNormal;	// Generic field outside valleys
  G4CMPEqEMField *EqValley1;		// Field for each of four Ge valleys
  G4CMPEqEMField *EqValley2;
  G4CMPEqEMField *EqValley3;
  G4CMPEqEMField *EqValley4;

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
};


#endif

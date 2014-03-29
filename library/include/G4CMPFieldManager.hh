// $Id$

#ifndef G4CMPFieldManager_h
#define G4CMPFieldManager_h 1

#include "globals.hh"
#include "G4FieldManager.hh"
#include "G4AffineTransform.hh"
#include <vector>

class G4CMPEqEMField;
class G4ChordFinder;
class G4ElectricField;
class G4EqMagElectricField;
class G4LatticeLogical;
class G4MagInt_Driver;
class G4MagIntegratorStepper;


class G4CMPFieldManager : public G4FieldManager {
public:
  G4CMPFieldManager(G4ElectricField *detectorField,
		    const G4LatticeLogical* lattice);

  ~G4CMPFieldManager();

  void ConfigureForTrack(const G4Track* aTrack);

private:
  const G4int stepperVars;
  const G4double stepperLength;

  // Non-const access to underlying field, base class doesn't provide
  G4ElectricField* myDetectorField;

  // NOTE: All pointers are kept in order to delete in dtor
  G4EqMagElectricField* EqNormal;	// Generic field outside valleys
  G4MagIntegratorStepper* normalStepper;
  G4MagInt_Driver* normalDriver;
  G4ChordFinder* normalChordFinder;

  const G4LatticeLogical* theLattice;
  size_t nValleys;
  std::vector<G4CMPEqEMField*> EqValley;	// Field for each Ge valley
  std::vector<G4MagIntegratorStepper*> valleyStepper;
  std::vector<G4MagInt_Driver*> valleyDriver;
  std::vector<G4ChordFinder*> valleyChordFinder;

  // Coordinate transformation from track, to pass to field and eq.-of-motion
  G4AffineTransform fLocalToGlobal;
};

#endif	/* G4CMPFieldManager_h */

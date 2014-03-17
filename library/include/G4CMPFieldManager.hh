// $Id$

#ifndef G4CMPFieldManager_h
#define G4CMPFieldManager_h 1

#include "globals.hh"
#include "G4FieldManager.hh"
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
};

#endif	/* G4CMPFieldManager_h */

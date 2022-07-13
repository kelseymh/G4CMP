/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
//
// 20140329  Pass G4CMP field pointer, which handles local/global transform
// 20140331  Use G4CMPEqEMField for everything, now handles holes; drop local
//	     caching of lattice.
// 20170525  Destructor should be virtual; add "rule of five" copy/move
// 20170801  Add counter to track instances of null-lattice, for reflections.
// 20210901  Add local verbosity flag for reporting diagnostics; use instead
//	     of G4CMP global setting.

#ifndef G4CMPFieldManager_h
#define G4CMPFieldManager_h 1

#include "globals.hh"
#include "G4FieldManager.hh"
#include <vector>

class G4CMPEqEMField;
class G4CMPLocalElectroMagField;
class G4ChordFinder;
class G4ElectroMagneticField;
class G4LatticeLogical;
class G4LatticePhysical;
class G4MagInt_Driver;
class G4MagIntegratorStepper;


class G4CMPFieldManager : public G4FieldManager {
public:
  G4CMPFieldManager(G4ElectroMagneticField *globalField, G4int vb=0);
  G4CMPFieldManager(G4CMPLocalElectroMagField *detectorField, G4int vb=0);
  virtual ~G4CMPFieldManager();

  // Base class does not allow copying
  G4CMPFieldManager(const G4CMPFieldManager&) = delete;
  G4CMPFieldManager(G4CMPFieldManager&&) = delete;
  G4CMPFieldManager& operator=(const G4CMPFieldManager&) = delete;
  G4CMPFieldManager& operator=(G4CMPFieldManager&&) = delete;

  // Run-time configuration
  void SetVerboseLevel(G4int vb) { verboseLevel = vb; }
  G4int GetVerboseLevel() const { return verboseLevel; }

  void ConfigureForTrack(const G4Track* aTrack);
  void SetChargeValleyForTrack(const G4LatticePhysical* lat, G4int valley);

private:
  G4int verboseLevel;		// For reporting diagnostic progress

  // Non-const access to underlying field, base class doesn't provide
  G4CMPLocalElectroMagField* myDetectorField;

  const G4int stepperVars;
  const G4double stepperLength;

  G4int latticeNulls;		// Count consective cases of no lattice
  const G4int maxLatticeNulls;	// Maximum allowed cases (reflection == 2)

  // NOTE: All pointers are kept in order to delete in dtor
  void CreateTransport();
  G4CMPEqEMField* theEqMotion;
  G4MagIntegratorStepper* theStepper;
  G4MagInt_Driver* theDriver;
  G4ChordFinder* theChordFinder;
};

#endif	/* G4CMPFieldManager_h */

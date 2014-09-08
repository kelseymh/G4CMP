#ifndef G4CMPConfigManager_hh
#define G4CMPConfigManager_hh 1

// $Id$
// File:  G4CMPConfigManager.hh
//
// Description:	Singleton container class for user configuration of G4CMP
//		applications at runtime.  Looks for environment variables
//		at initialization to set default values; active values may
//		be changed via macro commands (see G4CMPConfigMessenger).
//
// 20140904  Michael Kelsey

#include "globals.hh"

class G4CMPConfigMessenger;


class G4CMPConfigManager {
public:
  ~G4CMPConfigManager();	// Must be public for end-of-job cleanup

  // Access current values
  static G4int GetVerboseLevel()         { return Instance()->verbose; }
  static G4double GetVoltage()           { return Instance()->voltage; }
  static const G4String& GetEpotFile()   { return Instance()->Epot_file; }
  static const G4String& GetLatticeDir() { return Instance()->LatticeDir; }

  // Change values (e.g., via Messenger)
  static void SetVerboseLevel(G4int value) { Instance()->verbose = value; }
  static void SetVoltage(G4double value) { Instance()->voltage = value; }
  static void SetEpotFile(const G4String& name) { Instance()->Epot_file=name; }
  static void SetLatticeDir(const G4String& dir) { Instance()->LatticeDir=dir; }

private:
  G4CMPConfigManager();		// Singleton: only constructed on request

  static G4CMPConfigManager* Instance();   // Only needed by static accessors
  static G4CMPConfigManager* theInstance;

private:
  G4int verbose;		// Global verbosity (all processes, lattices)
  G4double voltage;		// Voltage setting for uniform electric field
  G4String Epot_file;		// Name of E-field file ($G4CMP_EPOT_FILE)
  G4String LatticeDir;		// Lattice data directory ($G4LATTICEDATA)

  G4CMPConfigMessenger* messenger;
};

#endif	/* G4CMPConfigManager_hh */

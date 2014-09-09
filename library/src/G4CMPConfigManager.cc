// $Id$
// File:  G4CMPConfigManager.cc
//
// Description:	Singleton container class for user configuration of G4CMP
//		applications at runtime.  Looks for environment variables
//		at initialization to set default values; active values may
//		be changed via macro commands (see G4CMPConfigMessenger).
//
// 20140904  Michael Kelsey

#include "G4CMPConfigManager.hh"
//*** #include "G4CMPConfigMessenger.hh"
#include <stdlib.h>


// Constructor and Singleton Initializer

G4CMPConfigManager* G4CMPConfigManager::theInstance = 0;

G4CMPConfigManager* G4CMPConfigManager::Instance() {
  if (!theInstance) theInstance = new G4CMPConfigManager;
  return theInstance;
}

G4CMPConfigManager::G4CMPConfigManager()
  : verbose(getenv("G4CMP_DEBUG")?atoi(getenv("G4CMP_DEBUG")):0),
    voltage(getenv("G4CMP_VOLTAGE")?strtod(getenv("G4CMP_VOLTAGE"),0):0.),
    Epot_file(getenv("G4CMP_EPOT_FILE")?getenv("G4CMP_EPOT_FILE"):"Epot_iZip4_small"),
    LatticeDir(getenv("G4LATTICEDATA")?getenv("G4LATTICEDATA"):"."),
    messenger(0 /* new G4CMPConfigMessenger(this) */) {;}

G4CMPConfigManager::~G4CMPConfigManager() {
  //*** delete messenger; messenger=0;
}

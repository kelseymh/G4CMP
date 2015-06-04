// $Id$
// File:  G4CMPConfigManager.cc
//
// Description:	Singleton container class for user configuration of G4CMP
//		applications at runtime.  Looks for environment variables
//		at initialization to set default values; active values may
//		be changed via macro commands (see G4CMPConfigMessenger).
//
// 20140904  Michael Kelsey
// 20141029  Force numerical voltage to correct units
// 20150603  Add parameter to limit reflections in DriftBoundaryProcess

#include "G4CMPConfigManager.hh"
#include "G4CMPConfigMessenger.hh"
#include "G4SystemOfUnits.hh"
#include <stdlib.h>


// Constructor and Singleton Initializer

G4CMPConfigManager* G4CMPConfigManager::theInstance = 0;

G4CMPConfigManager* G4CMPConfigManager::Instance() {
  if (!theInstance) theInstance = new G4CMPConfigManager;
  return theInstance;
}

G4CMPConfigManager::G4CMPConfigManager()
  : verbose(getenv("G4CMP_DEBUG")?atoi(getenv("G4CMP_DEBUG")):0),
    voltage(getenv("G4CMP_VOLTAGE")?strtod(getenv("G4CMP_VOLTAGE"),0)*volt:0.),
    stepScale(getenv("G4CMP_MIN_STEP")?strtod(getenv("G4CMP_MIN_STEP"),0):-1.),
    lukePhonons(getenv("G4CMP_LUKE_PHONONS")?strtod(getenv("G4CMP_LUKE_PHONONS"),0):0.),
    epotScale(getenv("G4CMP_EPOT_SCALE")?strtod(getenv("G4CMP_EPOT_SCALE"),0):1.),
    ehBounces(getenv("G4CMP_EH_BOUNCES")?atoi(getenv("G4CMP_EH_BOUNCES")):1),
    Epot_file(getenv("G4CMP_EPOT_FILE")?getenv("G4CMP_EPOT_FILE"):"Epot_iZip4_small"),
    LatticeDir(getenv("G4LATTICEDATA")?getenv("G4LATTICEDATA"):"."),
    Hit_file(getenv("G4CMP_HIT_FILE")?getenv("G4CMP_HIT_FILE"):"epositions.txt"),
    messenger(new G4CMPConfigMessenger(this)) {;}

G4CMPConfigManager::~G4CMPConfigManager() {
  delete messenger; messenger=0;
}

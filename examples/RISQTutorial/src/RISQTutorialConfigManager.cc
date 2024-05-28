/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
// File:  RISQTutorialConfigManager.cc
//
// Description:	Singleton container class for user configuration of G4CMP
//		phonon example. Looks for environment variables	at
//		initialization to set default values; active values may be
//		changed via macro commands (see RISQTutorialConfigMessenger).
//
// 20170816  M. Kelsey -- Extract hit filename from G4CMPConfigManager.

#include "RISQTutorialConfigManager.hh"
#include "RISQTutorialConfigMessenger.hh"
#include "G4RunManager.hh"
#include <stdlib.h>


// Constructor and Singleton Initializer

RISQTutorialConfigManager* RISQTutorialConfigManager::theInstance = 0;

RISQTutorialConfigManager* RISQTutorialConfigManager::Instance() {
  if (!theInstance) theInstance = new RISQTutorialConfigManager;
  return theInstance;
}

RISQTutorialConfigManager::RISQTutorialConfigManager()
  : Hit_file(getenv("G4CMP_HIT_FILE")?getenv("G4CMP_HIT_FILE"):"RISQTutorial_hits.txt"),
    Primary_file("RISQTutorial_primary.txt"),
    messenger(new RISQTutorialConfigMessenger(this)) {;}

RISQTutorialConfigManager::~RISQTutorialConfigManager() {
  delete messenger; messenger=0;
}


// Trigger rebuild of geometry if parameters change

void RISQTutorialConfigManager::UpdateGeometry() {
  G4RunManager::GetRunManager()->ReinitializeGeometry(true);
}

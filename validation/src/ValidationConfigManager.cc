/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file ValidationConfigManager.cc
/// \brief Description:	Singleton container class for user configuration of
//      G4CMP validation example. Looks for environment variables at
//		initialization to set default values; active values may be
//		changed via macro commands (see ValidationConfigMessenger).

#include "ValidationConfigManager.hh"
#include "ValidationConfigMessenger.hh"
#include "G4RunManager.hh"
#include <stdlib.h>


// Constructor and Singleton Initializer

ValidationConfigManager* ValidationConfigManager::theInstance = 0;

ValidationConfigManager* ValidationConfigManager::Instance() {
  if (!theInstance) theInstance = new ValidationConfigManager;
  return theInstance;
}

ValidationConfigManager::ValidationConfigManager()
  : Hit_file(getenv("G4CMP_HIT_FILE")?getenv("G4CMP_HIT_FILE"):"phonon_hits.txt"),
    stepFile("StepOutput.txt"),
    messenger(new ValidationConfigMessenger(this)) {;}

ValidationConfigManager::~ValidationConfigManager() {
  delete messenger; messenger=0;
}


// Trigger rebuild of geometry if parameters change

void ValidationConfigManager::UpdateGeometry() {
  G4RunManager::GetRunManager()->ReinitializeGeometry(true);
}

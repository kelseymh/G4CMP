/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// 20241024 Israel Hernandez -- IIT, QSC and Fermilab


#include "Caustic_PhononConfigManager.hh"
#include "Caustic_PhononConfigMessenger.hh"
#include "G4RunManager.hh"
#include <stdlib.h>


// Constructor and Singleton Initializer

Caustic_PhononConfigManager* Caustic_PhononConfigManager::theInstance = 0;

Caustic_PhononConfigManager* Caustic_PhononConfigManager::Instance() {
  if (!theInstance) theInstance = new Caustic_PhononConfigManager;
  return theInstance;
}

Caustic_PhononConfigManager::Caustic_PhononConfigManager()
  : Hit_file(getenv("G4CMP_HIT_FILE")?getenv("G4CMP_HIT_FILE"):"phonon_hits.txt"),
    messenger(new Caustic_PhononConfigMessenger(this))
{;}

Caustic_PhononConfigManager::~Caustic_PhononConfigManager() {
  delete messenger; messenger=0;
}


// Trigger rebuild of geometry if parameters change

void Caustic_PhononConfigManager::UpdateGeometry() {
  G4RunManager::GetRunManager()->ReinitializeGeometry(true);
}

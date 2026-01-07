/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
// File:  QuasiparticleConfigManager.cc
//
// Description:	Singleton container class for user configuration of G4CMP
//		quasiparticle example. Looks for environment variables	at
//		initialization to set default values; active values may be
//		changed via macro commands (see QuasiparticleConfigMessenger).
//
// 20170816  M. Kelsey -- Extract hit filename from G4CMPConfigManager.

#include "QuasiparticleConfigManager.hh"
#include "QuasiparticleConfigMessenger.hh"
#include "G4RunManager.hh"
#include <stdlib.h>


// Constructor and Singleton Initializer

QuasiparticleConfigManager* QuasiparticleConfigManager::theInstance = 0;

QuasiparticleConfigManager* QuasiparticleConfigManager::Instance() {
  if (!theInstance) theInstance = new QuasiparticleConfigManager;
  return theInstance;
}

QuasiparticleConfigManager::QuasiparticleConfigManager()
  : Hit_file(getenv("G4CMP_HIT_FILE")?getenv("G4CMP_HIT_FILE"):"phonon_hits.txt"),
    messenger(new QuasiparticleConfigMessenger(this)) {;}

QuasiparticleConfigManager::~QuasiparticleConfigManager() {
  delete messenger; messenger=0;
}


// Trigger rebuild of geometry if parameters change

void QuasiparticleConfigManager::UpdateGeometry() {
  G4RunManager::GetRunManager()->ReinitializeGeometry(true);
}

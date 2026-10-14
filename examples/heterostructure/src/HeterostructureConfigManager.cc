/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// /file  HeterostructureConfigManager.cc
/// /brief Implementation of the HeterostructureConfigManager class
///   Singleton container class for user configuration of G4CMP
///		Heterostructure example. Looks for environment variables	at
///		initialization to set default values; active values may be
///		changed via macro commands (see HeterostructureConfigMessenger).

// 20260815  Selby Q. Dang (Stanford/SLAC)

#include "HeterostructureConfigManager.hh"
#include "HeterostructureConfigMessenger.hh"
#include "G4RunManager.hh"
#include <stdlib.h>

// Constructor and Singleton Initializer

HeterostructureConfigManager* HeterostructureConfigManager::theInstance = 0;

HeterostructureConfigManager* HeterostructureConfigManager::Instance() {
  if (!theInstance) theInstance = new HeterostructureConfigManager;
  return theInstance;
}

HeterostructureConfigManager::HeterostructureConfigManager()
  : Hit_file(getenv("G4CMP_HIT_FILE")?getenv("G4CMP_HIT_FILE"):"Heterostructure_hits.txt"),
    fQAbsProbGeSi(0), fQAbsProbSiGe(0), fQReflProbGeSi(0), fQReflProbSiGe(0),
    messenger(new HeterostructureConfigMessenger(this)) {;}

HeterostructureConfigManager::~HeterostructureConfigManager() {
  delete messenger; messenger=0;
}


// Trigger rebuild of geometry if parameters change

void HeterostructureConfigManager::UpdateGeometry() {
  G4RunManager::GetRunManager()->ReinitializeGeometry(true);
}

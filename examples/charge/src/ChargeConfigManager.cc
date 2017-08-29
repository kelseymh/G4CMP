/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
// File:  ChargeConfigManager.cc
//
// Description:	Singleton container class for user configuration of G4CMP
//		phonon example. Looks for environment variables	at
//		initialization to set default values; active values may be
//		changed via macro commands (see ChargeConfigMessenger).
//
// 20170816  M. Kelsey -- Extract hit filename from G4CMPConfigManager.

#include "ChargeConfigManager.hh"
#include "ChargeConfigMessenger.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include <stdlib.h>


// Constructor and Singleton Initializer

ChargeConfigManager* ChargeConfigManager::theInstance = 0;

ChargeConfigManager* ChargeConfigManager::Instance() {
  if (!theInstance) theInstance = new ChargeConfigManager;
  return theInstance;
}

ChargeConfigManager::ChargeConfigManager()
  : voltage(getenv("G4CMP_VOLTAGE")?strtod(getenv("G4CMP_VOLTAGE"),0)*volt:0.),
    epotScale(getenv("G4CMP_EPOT_SCALE")?strtod(getenv("G4CMP_EPOT_SCALE"),0):1.),
    EPot_file(getenv("G4CMP_EPOT_FILE")?getenv("G4CMP_EPOT_FILE"):""),
    Hit_file(getenv("G4CMP_HIT_FILE")?getenv("G4CMP_HIT_FILE"):"charge_hits.txt"),
    millerH(0), millerK(0), millerL(0),
    messenger(new ChargeConfigMessenger(this)) {;}

ChargeConfigManager::~ChargeConfigManager() {
  delete messenger; messenger=0;
}


// Trigger rebuild of geometry if parameters change

void ChargeConfigManager::UpdateGeometry() {
  G4RunManager::GetRunManager()->ReinitializeGeometry(true);
}

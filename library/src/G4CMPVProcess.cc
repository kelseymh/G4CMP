/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPVProcess.cch
/// \brief Implementation of the G4CMPVProcess base class
//
// $Id$
//
// 20170601  New abstract base class for all G4CMP processes

#include "G4CMPVProcess.hh"
#include "G4CMPConfigManager.hh"


// Constructor and initialization only

G4CMPVProcess::G4CMPVProcess(const G4String& processName,
			     G4CMPProcessSubType stype)
  : G4VDiscreteProcess(processName, fPhonon), G4CMPProcessUtils() {
  verboseLevel = G4CMPConfigManager::GetVerboseLevel();
  SetProcessSubType(stype);
}

void G4CMPVProcess::StartTracking(G4Track* track) {
  G4VProcess::StartTracking(track);	// Apply base class actions
  LoadDataForTrack(track);
}

void G4CMPVProcess::EndTracking() {
  G4VProcess::EndTracking();		// Apply base class actions
  ReleaseTrack();
}

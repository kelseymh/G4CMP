/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef RISQTutorialConfigManager_hh
#define RISQTutorialConfigManager_hh 1

// $Id$
// File:  RISQTutorialConfigManager.hh
//
// Description:	Singleton container class for user configuration of G4CMP
//		phonon example. Looks for environment variables	at
//		initialization to set default values; active values may be
//		changed via macro commands (see RISQTutorialConfigMessenger).
//
// 20170816  M. Kelsey -- Extract hit filename from G4CMPConfigManager.

#include "globals.hh"

class RISQTutorialConfigMessenger;


class RISQTutorialConfigManager {
public:
  ~RISQTutorialConfigManager();	// Must be public for end-of-job cleanup
  static RISQTutorialConfigManager* Instance();   // Only needed by static accessors

  // Access current values
  static const G4String& GetHitOutput()  { return Instance()->Hit_file; }
  static const G4String& GetPrimaryOutput()  { return Instance()->Primary_file; }

  // Change values (e.g., via Messenger)
  static void SetHitOutput(const G4String& name)
    { Instance()->Hit_file=name; UpdateGeometry(); }

  // Change values (e.g., via Messenger)
  static void SetPrimaryOutput(const G4String& name)
    { Instance()->Hit_file=name; UpdateGeometry(); }

  
  static void UpdateGeometry();

private:
  RISQTutorialConfigManager();		// Singleton: only constructed on request
  RISQTutorialConfigManager(const RISQTutorialConfigManager&) = delete;
  RISQTutorialConfigManager(RISQTutorialConfigManager&&) = delete;
  RISQTutorialConfigManager& operator=(const RISQTutorialConfigManager&) = delete;
  RISQTutorialConfigManager& operator=(RISQTutorialConfigManager&&) = delete;

  static RISQTutorialConfigManager* theInstance;

private:
  G4String Hit_file;	// Output file of e/h hits ($G4CMP_HIT_FILE)
  G4String Primary_file;	// Output file of primaries

  RISQTutorialConfigMessenger* messenger;
};

#endif	/* RISQTutorialConfigManager_hh */

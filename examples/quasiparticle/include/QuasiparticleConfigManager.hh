/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef QuasiparticleConfigManager_hh
#define QuasiparticleConfigManager_hh 1

// $Id$
// File:  QuasiparticleConfigManager.hh
//
// Description:	Singleton container class for user configuration of G4CMP
//		quasiparticle example. Looks for environment variables	at
//		initialization to set default values; active values may be
//		changed via macro commands (see QuasiparticleConfigMessenger).
//
// 20170816  M. Kelsey -- Extract hit filename from G4CMPConfigManager.

#include "globals.hh"

class QuasiparticleConfigMessenger;


class QuasiparticleConfigManager {
public:
  ~QuasiparticleConfigManager();	// Must be public for end-of-job cleanup
  static QuasiparticleConfigManager* Instance();   // Only needed by static accessors

  // Access current values
  static const G4String& GetHitOutput()  { return Instance()->Hit_file; }

  // Change values (e.g., via Messenger)
  static void SetHitOutput(const G4String& name)
    { Instance()->Hit_file=name; UpdateGeometry(); }

  static void UpdateGeometry();

private:
  QuasiparticleConfigManager();		// Singleton: only constructed on request
  QuasiparticleConfigManager(const QuasiparticleConfigManager&) = delete;
  QuasiparticleConfigManager(QuasiparticleConfigManager&&) = delete;
  QuasiparticleConfigManager& operator=(const QuasiparticleConfigManager&) = delete;
  QuasiparticleConfigManager& operator=(QuasiparticleConfigManager&&) = delete;

  static QuasiparticleConfigManager* theInstance;

private:
  G4String Hit_file;	// Output file of e/h hits ($G4CMP_HIT_FILE)

  QuasiparticleConfigMessenger* messenger;
};

#endif	/* QuasiparticleConfigManager_hh */

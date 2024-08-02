/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef Caustic_PhononConfigManager_hh
#define Caustic_PhononConfigManager_hh 1
//20240110 Israel Hernandez -- Illinois Institute of Technology, Quantum Science Center and Fermilab

#include "globals.hh"

class Caustic_PhononConfigMessenger;


class Caustic_PhononConfigManager {
public:
  ~Caustic_PhononConfigManager();	// Must be public for end-of-job cleanup
  static Caustic_PhononConfigManager* Instance();   // Only needed by static accessors

  // Access current values
  static const G4String& GetHitOutput()  { return Instance()->Hit_file; }


  // Change values (e.g., via Messenger)
  static void SetHitOutput(const G4String& name)
    { Instance()->Hit_file=name; UpdateGeometry(); }


  static void UpdateGeometry();

private:
  Caustic_PhononConfigManager();		// Singleton: only constructed on request
  Caustic_PhononConfigManager(const Caustic_PhononConfigManager&) = delete;
  Caustic_PhononConfigManager(Caustic_PhononConfigManager&&) = delete;
  Caustic_PhononConfigManager& operator=(const Caustic_PhononConfigManager&) = delete;
  Caustic_PhononConfigManager& operator=(Caustic_PhononConfigManager&&) = delete;

  static Caustic_PhononConfigManager* theInstance;

private:
  G4String Hit_file;	// Output file 


  Caustic_PhononConfigMessenger* messenger;
};

#endif

/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef ChargeConfigManager_hh
#define ChargeConfigManager_hh 1

// $Id$
// File:  ChargeConfigManager.hh
//
// Description:	Singleton container class for user configuration of G4CMP
//		charge example. Looks for environment variables	at
//		initialization to set default values; active values may be
//		changed via macro commands (see ChargeConfigMessenger).
//
// 20170816  M. Kelsey -- Extract hit filename from G4CMPConfigManager.

#include "globals.hh"

class ChargeConfigMessenger;


class ChargeConfigManager {
public:
  ~ChargeConfigManager();	// Must be public for end-of-job cleanup
  static ChargeConfigManager* Instance();   // Only needed by static accessors

  // Access current values
  static G4double GetVoltage()           { return Instance()->voltage; }
  static G4double GetEPotScale()         { return Instance()->epotScale; }
  static const G4String& GetEPotFile()   { return Instance()->EPot_file; }
  static const G4String& GetHitOutput()  { return Instance()->Hit_file; }
  static G4int GetMillerH()		 { return Instance()->millerH; }
  static G4int GetMillerK()		 { return Instance()->millerK; }
  static G4int GetMillerL()		 { return Instance()->millerL; }

  static void GetMillerOrientation(G4int& h, G4int& k, G4int& l) {
    h = Instance()->millerH; k = Instance()->millerK; l = Instance()->millerL;
  }

  // Change values (e.g., via Messenger); require geometry to be rebuilt
  static void SetVoltage(G4double value)
    { Instance()->voltage = value; UpdateGeometry(); }
  static void SetEPotScale(G4double value)
    { Instance()->epotScale = value; UpdateGeometry(); }
  static void SetEPotFile(const G4String& name)
    { Instance()->EPot_file=name; UpdateGeometry(); }
  static void SetHitOutput(const G4String& name)
    { Instance()->Hit_file=name; UpdateGeometry(); }
  static void SetMillerOrientation(G4int h, G4int k, G4int l)
    { Instance()->millerH=h; Instance()->millerK=k, Instance()->millerL=l;
      UpdateGeometry();
    }

  static void UpdateGeometry();

private:
  ChargeConfigManager();		// Singleton: only constructed on request
  ChargeConfigManager(const ChargeConfigManager&) = delete;
  ChargeConfigManager(ChargeConfigManager&&) = delete;
  ChargeConfigManager& operator=(const ChargeConfigManager&) = delete;
  ChargeConfigManager& operator=(ChargeConfigManager&&) = delete;

  static ChargeConfigManager* theInstance;

private:
  G4double voltage;	// Uniform field voltage ($G4CMP_VOLTAGE)
  G4double epotScale;	// Scale factor for EPot ($G4CMP_EPOT_SCALE)
  G4String EPot_file;	// Name of E-field file ($G4CMP_EPOT_FILE)
  G4String Hit_file;	// Output file of e/h hits ($G4CMP_HIT_FILE)
  G4int millerH;	// Lattice orientation ($G4CMP_MILLER_H,_K,_L)
  G4int millerK;
  G4int millerL;

  ChargeConfigMessenger* messenger;
};

#endif	/* ChargeConfigManager_hh */

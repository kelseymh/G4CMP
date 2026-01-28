/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file:  ValidationConfigManager.hh
/// \brief Description:	Singleton container class for user configuration of
///     G4CMP validation example. Looks for environment variables at
///		initialization to set default values; active values may be
///		changed via macro commands (see ValidationConfigMessenger).
//
//  20260109  M. Kelsey -- Missing Geant4 type declarations

#ifndef ValidationConfigManager_hh
#define ValidationConfigManager_hh 1

#include "G4Types.hh"
#include "G4String.hh"

class ValidationConfigMessenger;


class ValidationConfigManager {
public:
  ~ValidationConfigManager();	// Must be public for end-of-job cleanup
  static ValidationConfigManager* Instance();   // Only needed by static accessors

  // Access current values
  static const G4String& GetHitOutput()  { return Instance()->Hit_file; }
  static const G4int& GetGeometryID() { return Instance()->geometryID; }
  static const G4String& GetStepFileName() { return Instance()->stepFile; }
  
  // Change values (e.g., via Messenger)
  static void SetHitOutput(const G4String& name)
  { Instance()->Hit_file=name; UpdateGeometry(); }
  static void SetGeometryID(const G4int& ID)
  { Instance()->geometryID=ID; UpdateGeometry(); } 
  static void SetStepFile(const G4String& name)
  { Instance()->stepFile=name; }
  
  static void UpdateGeometry();

private:
  ValidationConfigManager();		// Singleton: only constructed on request
  ValidationConfigManager(const ValidationConfigManager&) = delete;
  ValidationConfigManager(ValidationConfigManager&&) = delete;
  ValidationConfigManager& operator=(const ValidationConfigManager&) = delete;
  ValidationConfigManager& operator=(ValidationConfigManager&&) = delete;

  static ValidationConfigManager* theInstance;

private:
  G4String Hit_file;	// Output file of e/h hits ($G4CMP_HIT_FILE)
  G4int geometryID;
  G4String stepFile;

  ValidationConfigMessenger* messenger;
};

#endif	/* ValidationConfigManager_hh */

/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPVScatteringRate.hh
/// \brief Definition of the G4CMPVScatteringRate base class.  This class
///	   provides an interface to implement calculations of scattering
///	   rate (either phenomenological or theoretical) for phonons or
///	   charge carriers.
//
// $Id$

#ifndef G4CMPVScatteringRate_hh
#define G4CMPVScatteringRate_hh 1

#include "globals.hh"
#include "G4CMPConfigManager.hh"

class G4Track;


class G4CMPVScatteringRate {
public:
  G4CMPVScatteringRate(const G4String& theName)
    : name(theName), verboseLevel(G4CMPConfigManager::GetVerboseLevel()) {;}

  virtual ~G4CMPVScatteringRate() {;}

  // Get scattering rate from track kinematics; subclasses MUST IMPLEMENT
  virtual G4double Rate(const G4Track& aTrack) const = 0;

  // Additional interface which call back to above; should not be overridden
  G4double Rate(const G4Track* aTrack) const { return Rate(*aTrack); }

  // General configuration
  void SetVerboseLevel(G4int vb) { verboseLevel = vb; }
  G4int GetVerboseLevel() const { return verboseLevel; }

  const G4String& GetName() const { return name; }

protected:
  G4String name;		// Accessible for use by subclasses
  G4int verboseLevel;
};

#endif	/* G4CMPVScatteringRate_hh */

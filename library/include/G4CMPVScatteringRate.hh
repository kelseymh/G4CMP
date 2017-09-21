/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
//
/// \file library/include/G4CMPVScatteringRate.hh
/// \brief Definition of the G4CMPVScatteringRate base class.  This class
///	   provides an interface to implement calculations of scattering
///	   rate (either phenomenological or theoretical) for phonons or
///	   charge carriers.  Subclasses may specified ctor argument if
///	   process should be forced.
//
// 20170815  Inherit from G4CMPProcessUtils here, instead of in subclasses
// 20170919  Add "threshold finder" interface, for use with IV and Luke

#ifndef G4CMPVScatteringRate_hh
#define G4CMPVScatteringRate_hh 1

#include "globals.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPProcessUtils.hh"

class G4Track;


class G4CMPVScatteringRate : public G4CMPProcessUtils {
public:
  G4CMPVScatteringRate(const G4String& theName, G4bool force=false)
    : G4CMPProcessUtils(),
      verboseLevel(G4CMPConfigManager::GetVerboseLevel()),
      name(theName), isForced(force) {;}

  virtual ~G4CMPVScatteringRate() {;}

  // Get scattering rate from track kinematics; subclasses MUST IMPLEMENT
  virtual G4double Rate(const G4Track& aTrack) const = 0;

  // Additional interface which call back to above; should not be overridden
  G4double Rate(const G4Track* aTrack) const { return Rate(*aTrack); }

  // Interface to identify energy thresholds (for IV, Luke subclasses)
  virtual G4double Threshold(G4double /*Eabove*/=0.) const { return 0.; }

  // Flag if interaction should be forced (subclasses should set flag)

  G4bool IsForced() { return isForced; }

  // General configuration
  void SetVerboseLevel(G4int vb) { verboseLevel = vb; }
  G4int GetVerboseLevel() const { return verboseLevel; }

  const G4String& GetName() const { return name; }

protected:
  G4int verboseLevel;		// Accessible for use by subclasses
  G4String name;		// For diagnostic output if desired
  G4bool isForced;		// Flag 'true' if process should be forced
};

#endif	/* G4CMPVScatteringRate_hh */

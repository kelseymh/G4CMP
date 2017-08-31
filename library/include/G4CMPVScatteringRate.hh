/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

<<<<<<< HEAD
// $Id$
//
=======
>>>>>>> G4CMP-18
/// \file library/include/G4CMPVScatteringRate.hh
/// \brief Definition of the G4CMPVScatteringRate base class.  This class
///	   provides an interface to implement calculations of scattering
///	   rate (either phenomenological or theoretical) for phonons or
///	   charge carriers.  Subclasses may specified ctor argument if
///	   process should be forced.
//
<<<<<<< HEAD
// 20170815  Inherit from G4CMPProcessUtils here, instead of in subclasses
=======
// $Id$
>>>>>>> G4CMP-18

#ifndef G4CMPVScatteringRate_hh
#define G4CMPVScatteringRate_hh 1

#include "globals.hh"
#include "G4CMPConfigManager.hh"
<<<<<<< HEAD
#include "G4CMPProcessUtils.hh"
=======
>>>>>>> G4CMP-18

class G4Track;


<<<<<<< HEAD
class G4CMPVScatteringRate : public G4CMPProcessUtils {
public:
  G4CMPVScatteringRate(const G4String& theName, G4bool force=false)
    : G4CMPProcessUtils(),
      verboseLevel(G4CMPConfigManager::GetVerboseLevel()),
=======
class G4CMPVScatteringRate {
public:
  G4CMPVScatteringRate(const G4String& theName, G4bool force=false)
    : verboseLevel(G4CMPConfigManager::GetVerboseLevel()),
>>>>>>> G4CMP-18
      name(theName), isForced(force) {;}

  virtual ~G4CMPVScatteringRate() {;}

  // Get scattering rate from track kinematics; subclasses MUST IMPLEMENT
  virtual G4double Rate(const G4Track& aTrack) const = 0;

  // Additional interface which call back to above; should not be overridden
  G4double Rate(const G4Track* aTrack) const { return Rate(*aTrack); }

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

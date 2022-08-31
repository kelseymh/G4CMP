/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
//
/// \file  library/include/G4CMPVElectrodePattern.hh
/// \brief Abstract base class to define complex electrode layouts
//
// 20160831  M. Kelsey -- Add optional electrode geometry class
// 20160904  M. Kelsey -- Pass ref to concrete G4ParticleChange
// 20160906  M. Kelsey -- Add function handle constness of material table
// 20170525  M. Kelsey -- Add "rule of five" default copy/move operators
// 20170627  M. Kelsey -- Inherit from G4CMPProcessUtils
// 20200601  G4CMP-207: Require Clone() functions from sublcasses for copying

#ifndef G4CMPVElectrodePattern_h
#define G4CMPVElectrodePattern_h 1

#include "globals.hh"
#include "G4CMPProcessUtils.hh"
#include "G4MaterialPropertiesTable.hh"

class G4CMPSurfaceProperty;
class G4ParticleChange;
class G4Step;
class G4Track;


class G4CMPVElectrodePattern : public G4CMPProcessUtils {
public:
  G4CMPVElectrodePattern() : verboseLevel(0) {;}
  virtual ~G4CMPVElectrodePattern() {;}

  // Use default copy/move operators
  G4CMPVElectrodePattern(const G4CMPVElectrodePattern&) = default;
  G4CMPVElectrodePattern(G4CMPVElectrodePattern&&) = default;
  G4CMPVElectrodePattern& operator=(const G4CMPVElectrodePattern&) = default;
  G4CMPVElectrodePattern& operator=(G4CMPVElectrodePattern&&) = default;

  // Subclasses MUST implement this to make thread-local copies
  // NOTE: This default is NOT THREAD SAFE; it reproduces old behaviour
  virtual G4CMPVElectrodePattern* Clone() const {
    return const_cast<G4CMPVElectrodePattern*>(this);
  }

  // Subclasses may use verbosity level for diagnostics
  void SetVerboseLevel(G4int vb) { verboseLevel = vb; }

  // Local copy of properties stored automatically by G4CMPSurfaceProperty
  void UseSurfaceTable(G4MaterialPropertiesTable* surfProp) {
    theSurfaceTable = surfProp;
  }

  // Subclass MUST implement this to return true/false depending on position
  virtual G4bool IsNearElectrode(const G4Step& aStep) const = 0;

  // Subclass MAY implement this to deposit energy from track into electrode
  // Return is not necessary: aParticleChange may be altered in situ
  virtual void AbsorbAtElectrode(const G4Track&, const G4Step&,
                                 G4ParticleChange&) const {;}

protected:
  // Handles casting table to non-const for access
  G4double GetMaterialProperty(const G4String& key) const;

  G4int verboseLevel;
  G4MaterialPropertiesTable* theSurfaceTable;	// Not owned, can't be const
};

#endif	/* G4CMPVElectrodePattern_h */

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

#ifndef G4CMPVElectrodePattern_h
#define G4CMPVElectrodePattern_h 1

class G4CMPSurfaceProperty;
class G4Step;
class G4Track;
class G4VParticleChange;


class G4CMPVElectrodePattern {
public:
  G4CMPVElectrodePattern() : verboseLevel(0) {;}
  virtual ~G4CMPVElectrodePattern() {;}

  // Subclasses may use verbosity level for diagnostics
  void SetVerboseLevel(G4int vb) { verboseLevel = vb; }

  void UseSurfaceTable(const G4MaterialPropertiesTable& surfProp) {
    theSurfaceTable = surfProp;
  }

  // Subclass MUST implement this to return true/false depending on position
  virtual G4bool IsNearElectrode(const G4Step& aStep) const = 0;

  // Subclass MAY implement this to deposit energy from track into electrode
  // Return is not necessary: aParticleChange may be altered in situ
  virtual void AbsorbAtElectrode(const G4Track& aTrack, const G4Step& aStep,
				 G4VParticleChange* aParticleChange) const {;}

protected:
  G4int verboseLevel;
  G4MaterialPropertiesTable theSurfaceTable;
};

#endif	/* G4CMPVElectrodePattern_h */

/***********************************************************************\
 *  * This software is licensed under the terms of the GNU General Public *
 *   * License version 3 or later. See G4CMP/LICENSE for the full license. *
 *   \***********************************************************************/

/// \file library/include/G4CMPKaplanUtils.hh
//
//


#ifndef G4CMPKaplanUtils_hh
#define G4CMPKaplanUtils_hh 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include <utility>

class G4CMPProcessUtils;
class G4CMPSurfaceProperty;
class G4MaterialPropertiesTable;
class G4ParticleChange;
class G4Step;
class G4Track;
class G4VPhysicalVolume;
class G4VProcess;

class G4CMPKaplanUtils {
public:

  G4CMPKaplanUtils(){;}
  virtual ~G4CMPKaplanUtils(){;}

  G4CMPKaplanUtils(const G4CMPKaplanUtils&) = default;
  G4CMPKaplanUtils(G4CMPKaplanUtils&&) = default;
  G4CMPKaplanUtils& operator=(const G4CMPKaplanUtils&) = default;
  G4CMPKaplanUtils& operator=(G4CMPKaplanUtils&&) = default;
 
  virtual void SetVerboseLevel(G4int vb) { buVerboseLevel = vb; }


  // Compute quasiparticle energy distribution from broken Cooper pair.
  G4double QPEnergyRand(G4double Energy) const;
  G4double QPEnergyPDF(G4double E, G4double x) const;

  // Compute phonon energy distribution from quasiparticle in superconductor.
  G4double PhononEnergyRand(G4double Energy) const;
  G4double PhononEnergyPDF(G4double E, G4double x) const;

  G4bool IsSubgap(G4double energy) const { return energy < 2.*gapEnergy; }
  G4double getGapEnergy() const { return gapEnergy; }

private:
  G4int buVerboseLevel;			// For local use; name avoids collisions
  G4double gapEnergy;

};

#endif

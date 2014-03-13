//
/// \file library/include/G4CMPProcessUtils.hh
/// \brief Definition of the G4CMPProcessUtils class
///   Provides useful general functions to fetch and store lattice, access
///   and apply lattice parameters for phonons and charge carriers.
///
///   Use via multiple inheritance with concrete or base process classes
//
// $Id$
//

#ifndef G4CMPProcessUtils_hh
#define G4CMPProcessUtils_hh 1

#include "globals.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

class G4LatticePhysical;
class G4Track;
class G4PhononTrackMap;
class G4CMPValleyTrackMap;


class G4CMPProcessUtils {
public:
  G4CMPProcessUtils();
  ~G4CMPProcessUtils();

  // Fetch lattice (if any) for volume occupied by track
  void LoadDataForTrack(const G4Track* track);
  void ReleaseTrack();

  // Map phonon types to polarization index
  G4int GetPolarization(const G4Track& track) const;
  G4int GetPolarization(const G4Track* track) const {
    return GetPolarization(*track);
  }

  // Generate random phonon polarization from density of states
  // Values passed may be zero to suppress particular states
  G4int ChoosePolarization(G4double Ldos, G4double STdos, G4double FTdos) const;

  // Map charge carrier momentum to valley index
  G4int GetValleyIndex(const G4Track& track) const;
  G4int GetValleyIndex(const G4Track* track) const {
    return GetValleyIndex(*track);
  }

  const G4RotationMatrix& GetValley(const G4Track& track) const;
  const G4RotationMatrix& GetValley(const G4Track* track) const {
    return GetValley(*track); 
  }

  // Generate random valley for charge carrier
  G4int ChooseValley() const;

  // Construct new phonon track with correct momentum, position, etc.
  G4Track* CreatePhonon(G4int polarization, const G4ThreeVector& K,
			G4double energy) const;

  // Construct new electron or hole track with correct conditions
  G4Track* CreateChargeCarrier(G4int charge, G4int valley,
			       const G4ThreeVector& p, G4double energy) const;

protected:
  const G4LatticePhysical* theLattice;

  G4PhononTrackMap* trackKmap;		// For convenient access by processes
  G4CMPValleyTrackMap* trackVmap;

private:
  const G4Track* currentTrack;		// For use by Start/EndTracking

  // hide assignment operators as private 
  G4CMPProcessUtils(G4CMPProcessUtils&);
  G4CMPProcessUtils& operator=(const G4CMPProcessUtils& right);
};

#endif	/* G4CMPProcessUtils_hh */

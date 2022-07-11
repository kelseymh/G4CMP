/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPProcessUtils.hh
/// \brief Definition of the G4CMPProcessUtils class
///   Provides useful general functions to fetch and store lattice, access
///   and apply lattice parameters for phonons and charge carriers.
///
///   Use via multiple inheritance with concrete or base process classes
//
// $Id$
//
// 20140321  Move lattice-based placement transformations here, via Touchable
// 20140407  Add functions for phonon generation in Luke scattering
// 20140412  Add manual configuration options
// 20140509  Add ChoosePolarization() which uses DOS values from lattice
// 20150112  Add GetCurrentValley() function to get valley of current track,
//	     protect valley functions against null pointers
// 20150309  Add Create*() functions which take position and energy arguments
//	     (for use with AlongStepDoIt() actions).
// 20151209  Replace trackMap classes with an AuxiliaryInformation class.
// 20160107  Add GetVelocity functions.
// 20160610  Add accessor for auxiliary information
// 20160625  Add accessors for particle identification
// 20160825  Add assignment operators for cross-process configuration;
//	     move track identification functions to G4CMPUtils
// 20160906  Make GetSurfaceNormal() const.
// 20161004  Add new ChangeValley() function to avoid null selection
// 20170525  Drop explicit copy constructors; let compiler do the work
// 20170602  Local track identification functions apply to current track only
// 20170603  Drop deprecated functions; don't deprecate transforms.
// 20170620  Drop local caching of transforms; call through to G4CMPUtils.
// 20170806  Move ChargeCarrierTimeStep() here from DriftProcess.
// 20201111  Add MakePhononEnergy() which takes wave vector directly
// 20201124  Change argument name in MakeGlobalRecoil() to 'krecoil' (track)
// 20201223  Add FindNearestValley() function to align electron momentum.

#ifndef G4CMPProcessUtils_hh
#define G4CMPProcessUtils_hh 1

#include "globals.hh"
#include "G4AffineTransform.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"

class G4CMPDriftTrackInfo;
class G4CMPPhononTrackInfo;
class G4CMPVTrackInfo;
class G4LatticePhysical;
class G4ParticleDefinition;
class G4VPhysicalVolume;
class G4VTouchable;


class G4CMPProcessUtils {
public:
  G4CMPProcessUtils();
  virtual ~G4CMPProcessUtils();

  // Assignment operators allow dependent configuration
  G4CMPProcessUtils(const G4CMPProcessUtils&) = default;
  G4CMPProcessUtils(G4CMPProcessUtils&&) = default;
  G4CMPProcessUtils& operator=(const G4CMPProcessUtils&) = default;
  G4CMPProcessUtils& operator=(G4CMPProcessUtils&&) = default;

  // Configure for current track
  virtual void LoadDataForTrack(const G4Track* track);
  virtual void SetCurrentTrack(const G4Track* track);
  virtual void SetLattice(const G4Track* track);

  virtual void ReleaseTrack();
  // NOTE:  Subclasses may overload these, but be sure to callback to base

  // Identify current track type to simplify some conditionals
  G4bool IsPhonon() const;
  G4bool IsElectron() const;
  G4bool IsHole() const;
  G4bool IsChargeCarrier() const;

  // Set configuration manually, without a track
  virtual void FindLattice(const G4VPhysicalVolume* volume);
  virtual void SetLattice(const G4LatticePhysical* lat) { theLattice = lat; }
  virtual const G4LatticePhysical* GetLattice() const { return theLattice; }

  // Convert global to local coordinates with respect to current track
  G4ThreeVector GetLocalDirection(const G4ThreeVector& dir) const;
  G4ThreeVector GetLocalPosition(const G4ThreeVector& pos) const;
  void RotateToLocalDirection(G4ThreeVector& dir) const;
  void RotateToLocalPosition(G4ThreeVector& pos) const;

  // Convert local to global coordinates with respect to current track
  G4ThreeVector GetGlobalDirection(const G4ThreeVector& dir) const;
  G4ThreeVector GetGlobalPosition(const G4ThreeVector& pos) const;
  void RotateToGlobalDirection(G4ThreeVector& dir) const;
  void RotateToGlobalPosition(G4ThreeVector& pos) const;

  // Convenience functions to get local position, momentum, velocity from track
  G4ThreeVector GetLocalPosition(const G4Track& track) const;
  G4ThreeVector GetLocalPosition(const G4Track* track) const {
    return GetLocalPosition(*track);
  }

  void GetLocalPosition(const G4Track& track, G4double pos[3]) const;
  void GetLocalPosition(const G4Track* track, G4double pos[3]) const {
    return GetLocalPosition(*track, pos);
  }

  G4ThreeVector GetLocalMomentum(const G4Track& track) const;
  G4ThreeVector GetLocalMomentum(const G4Track* track) const {
    return GetLocalMomentum(*track);
  }

  void GetLocalMomentum(const G4Track& track, G4double mom[3]) const;
  void GetLocalMomentum(const G4Track* track, G4double mom[3]) const {
    return GetLocalMomentum(*track, mom);
  }

  G4ThreeVector GetLocalVelocityVector(const G4Track& track) const;
  G4ThreeVector GetLocalVelocityVector(const G4Track* track) const {
    return GetLocalVelocityVector(*track);
  }

  void GetLocalVelocityVector(const G4Track& track, G4double vel[3]) const;
  void GetLocalVelocityVector(const G4Track* track, G4double vel[3]) const {
    return GetLocalVelocityVector(*track, vel);
  }

  G4ThreeVector GetLocalWaveVector(const G4Track& track) const;
  G4ThreeVector GetLocalWaveVector(const G4Track* track) const {
    return GetLocalWaveVector(*track);
  }

  // Convenience functions to get position, momentum, velocity from track
  G4ThreeVector GetGlobalPosition(const G4Track& track) const;
  G4ThreeVector GetGlobalPosition(const G4Track* track) const {
    return GetGlobalPosition(*track);
  }

  void GetGlobalPosition(const G4Track& track, G4double pos[3]) const;
  void GetGlobalPosition(const G4Track* track, G4double pos[3]) const {
    return GetGlobalPosition(*track, pos);
  }

  G4ThreeVector GetGlobalMomentum(const G4Track& track) const;
  G4ThreeVector GetGlobalMomentum(const G4Track* track) const {
    return GetGlobalMomentum(*track);
  }

  void GetGlobalMomentum(const G4Track& track, G4double mom[3]) const;
  void GetGlobalMomentum(const G4Track* track, G4double mom[3]) const {
    return GetGlobalMomentum(*track, mom);
  }

  // These functions only exist to make a uniform interface for getting
  // track info in a G4CMP physics process class.
  G4double CalculateVelocity(const G4Track& track) const {
    return track.CalculateVelocity();
  }

  G4double CalculateVelocity(const G4Track* track) const {
    return CalculateVelocity(*track);
  }

  G4double GetVelocity(const G4Track& track) const {
    return CalculateVelocity(track);
  }

  G4double GetVelocity(const G4Track* track) const {
    return CalculateVelocity(*track);
  }

  G4ThreeVector GetGlobalVelocityVector(const G4Track& track) const;
  G4ThreeVector GetGlobalVelocityVector(const G4Track* track) const {
    return GetGlobalVelocityVector(*track);
  }

  void GetGlobalVelocityVector(const G4Track& track, G4double vel[3]) const;
  void GetGlobalVelocityVector(const G4Track* track, G4double vel[3]) const {
    return GetGlobalVelocityVector(*track, vel);
  }

  // Convenience functions to get position, momentum, velocity from track

  G4double GetKineticEnergy(const G4Track& track) const;
  G4double GetKineticEnergy(const G4Track* track) const {
    return GetKineticEnergy(*track);
  }

  // Map phonon types to polarization index
  G4int GetPolarization(const G4Track& track) const;
  inline G4int GetPolarization(const G4Track* track) const {
    return GetPolarization(*track);
  }

  // Generate random phonon polarization from density of states
  // Values passed may be zero to suppress particular states
  G4int ChoosePhononPolarization() const;		// Use DOS values from lattice

  // Map charge carrier momentum to valley index
  G4int GetValleyIndex(const G4Track& track) const;
  inline G4int GetValleyIndex(const G4Track* track) const {
    return (track ? GetValleyIndex(*track) : -1);
  }

  const G4RotationMatrix& GetValley(const G4Track& track) const;
  inline const G4RotationMatrix& GetValley(const G4Track* track) const {
    return (track ? GetValley(*track) : G4RotationMatrix::IDENTITY); 
  }

  // Generate random valley for charge carrier
  G4int ChangeValley(G4int valley) const;	// Excludes input valley

  // Find valley which aligns most closely with _local_ direction vector
  // NOTE: Passed by value to allow for internal manipulation
  G4int FindNearestValley(G4ThreeVector dir) const;
  G4int FindNearestValley(const G4Track& track) const;
  G4int FindNearestValley(const G4Track* track) const {
    return (track ? FindNearestValley(*track) : -1);
  }
  G4int FindNearestValley() const {		// Use current track
    return FindNearestValley(currentTrack);
  }

  // Generate direction angle for phonon generated in Luke scattering
  G4double MakePhononTheta(G4double k, G4double ks) const;
  G4double MakePhononEnergy(G4double k, G4double ks, G4double th_phonon) const;
  G4double MakePhononEnergy(G4double q) const;

  // Compute direction angle for recoiling charge carrier
  G4double MakeRecoilTheta(G4double k, G4double ks, G4double th_phonon) const;

  // Convert K_HV wave vector to track momentum
  void MakeGlobalRecoil(G4ThreeVector& krecoil) const;

  // Compute time between scatters/emissions for moving charge carrier
  // Parameters are "Mach number" (ratio to sound speed) and scattering length
  G4double ChargeCarrierTimeStep(G4double mach, G4double l0) const;

protected:
  const G4LatticePhysical* theLattice;	// For convenient access by processes

  const G4Track* GetCurrentTrack() const { return currentTrack; }
  const G4VPhysicalVolume* GetCurrentVolume() const { return currentVolume; }
  const G4ParticleDefinition* GetCurrentParticle() const;
  const G4VTouchable* GetCurrentTouchable() const;

  G4int GetCurrentValley() const { return GetValleyIndex(currentTrack); }

private:
  const G4Track* currentTrack;		// For use by Start/EndTracking
  const G4VPhysicalVolume* currentVolume;
};

#endif	/* G4CMPProcessUtils_hh */

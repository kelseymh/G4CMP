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

#ifndef G4CMPProcessUtils_hh
#define G4CMPProcessUtils_hh 1

#include "globals.hh"
#include "G4AffineTransform.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"

class G4CMPTrackInformation;
class G4LatticePhysical;
class G4ParticleDefinition;
class G4VPhysicalVolume;
class G4VTouchable;


class G4CMPProcessUtils {
public:
  G4CMPProcessUtils();
  virtual ~G4CMPProcessUtils();

  // Assignment operators allow dependent configuration
  G4CMPProcessUtils(G4CMPProcessUtils&);
  G4CMPProcessUtils& operator=(const G4CMPProcessUtils& right);

  // Configure for current track
  virtual void LoadDataForTrack(const G4Track* track);
  virtual void ReleaseTrack();
  // NOTE:  Subclasses may overload these, but be sure to callback to base

  // Identify track type to simplify some conditionals
  G4bool IsPhonon(const G4Track* track) const;
  G4bool IsElectron(const G4Track* track) const;
  G4bool IsHole(const G4Track* track) const;
  G4bool IsChargeCarrier(const G4Track* track) const;

  // Set configuration manually, without a track
  virtual void FindLattice(const G4VPhysicalVolume* volume);
  virtual void SetLattice(const G4LatticePhysical* lat) { theLattice = lat; }
  virtual const G4LatticePhysical* GetLattice() const { return theLattice; }

  virtual void SetTransforms(const G4VTouchable* touchable);
  virtual void SetTransforms(const G4RotationMatrix* rot,
			     const G4ThreeVector& trans);

  // Attach or retrieve auxiliary information object for track
  // NOTE:  This will cast the track to non-const if required
  G4CMPTrackInformation* AttachTrackInfo(const G4Track* track) const;

  // Extract auxiliary information for track (current track if none given)
  G4CMPTrackInformation* GetTrackInfo(const G4Track* track=0) const;
  G4CMPTrackInformation* GetTrackInfo(const G4Track& track) const {
    return GetTrackInfo(&track);
  }

  // Convert global to local coordinates
  inline G4ThreeVector GetLocalDirection(const G4ThreeVector& dir) const {
    return fGlobalToLocal.TransformAxis(dir);
  }

  inline G4ThreeVector GetLocalPosition(const G4ThreeVector& pos) const {
    return fGlobalToLocal.TransformPoint(pos);
  }

  inline void RotateToLocalDirection(G4ThreeVector& dir) const {
    fGlobalToLocal.ApplyAxisTransform(dir);
  }

  inline void RotateToLocalPosition(G4ThreeVector& pos) const {
    fGlobalToLocal.ApplyPointTransform(pos);
  }

  // Convert local to global coordinates
  inline G4ThreeVector GetGlobalDirection(const G4ThreeVector& dir) const {
    return fLocalToGlobal.TransformAxis(dir);
  }

  inline G4ThreeVector GetGlobalPosition(const G4ThreeVector& pos) const {
    return fLocalToGlobal.TransformPoint(pos);
  }

  inline void RotateToGlobalDirection(G4ThreeVector& dir) const {
    fLocalToGlobal.ApplyAxisTransform(dir);
  }

  inline void RotateToGlobalPosition(G4ThreeVector& pos) const {
    fLocalToGlobal.ApplyPointTransform(pos);
  }

  // Used especially in boundary processes
  G4ThreeVector GetSurfaceNormal(const G4Step& aStep) const;

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
  G4int ChooseValley() const;
  G4int ChangeValley(G4int valley) const;	// Excludes input valley

  // Generate direction angle for phonon generated in Luke scattering
  G4double MakePhononTheta(G4double k, G4double ks) const;
  G4double MakePhononEnergy(G4double k, G4double ks, G4double th_phonon) const;

  // Compute direction angle for recoiling charge carrier
  G4double MakeRecoilTheta(G4double k, G4double ks, G4double th_phonon) const;

  void MakeLocalPhononK(G4ThreeVector& kphonon) const;
  void MakeGlobalPhononK(G4ThreeVector& kphonon) const;

  void MakeGlobalRecoil(G4ThreeVector& kphonon) const;

  // Compute a Lambertian distribution for reflected phonons
  G4ThreeVector LambertReflection(const G4ThreeVector& surfNorm);

  // Model Kaplan phonon-quasiparticle interactions in superconductor sensors
  G4double KaplanPhononQP(G4double energy, G4MaterialPropertiesTable* prop,
                          std::vector<G4double>& reflectedEnergies);

  // Compute the probability of a phonon reentering the crystal
  G4double CalcEscapeProbability(G4double energy, G4double thicknessFrac,
                                 G4MaterialPropertiesTable* prop);

  // Model the phonons breaking Cooper pairs into quasiparticles
  G4double CalcQPEnergies(G4double gapEnergy, G4double lowQPLimit,
                          std::vector<G4double>& phonEnergies,
                          std::vector<G4double>& qpEnergies);

  // Model the quasiparticles emitting phonons in the superconductor
  G4double CalcPhononEnergies(G4double gapEnergy, G4double lowQPLimit,
                              std::vector<G4double>& phonEnergies,
                              std::vector<G4double>& qpEnergies);

  // Calculate energies of phonon tracks that have reentered the crystal
  void CalcReflectedPhononEnergies(G4MaterialPropertiesTable* prop,
                                   std::vector<G4double>& phonEnergies,
                                   std::vector<G4double>& reflectedEnergies);

  // Compute quasiparticle energy distribution from broken Cooper pair
  G4double QPEnergyRand(G4double gapEnergy, G4double Energy);

  // Compute phonon energy distribution from quasiparticle in superconductor
  G4double PhononEnergyRand(G4double gapEnergy, G4double& Energy);


  // Construct new phonon or charge carrier track
  G4Track* CreateTrack(G4ParticleDefinition* pd, const G4ThreeVector& K,
		       G4double energy) const;

  G4Track* CreateTrack(G4ParticleDefinition* pd, const G4ThreeVector& K,
		       G4double energy, const G4ThreeVector& pos) const;

  // Construct new phonon track with correct momentum, position, etc.
  G4Track* CreatePhonon(G4int polarization, const G4ThreeVector& K,
			G4double energy) const;

  G4Track* CreatePhonon(G4int polarization, const G4ThreeVector& K,
            G4double energy, const G4ThreeVector& pos) const;

  G4Track* CreatePhononInFromBoundary(G4int polarization,
                                      const G4ThreeVector& K,
                                      G4double energy) const;

  // Construct new electron or hole track with correct conditions
  G4Track* CreateChargeCarrier(G4int charge, G4int valley,
			       const G4ThreeVector& p) const;

  G4Track* CreateChargeCarrier(G4int charge, G4int valley, G4double Ekin,
			       const G4ThreeVector& pdir,
			       const G4ThreeVector& pos) const;

  G4Track* CreateChargeCarrier(G4int charge, G4int valley,
			       const G4ThreeVector& p,
                   const G4ThreeVector& pos) const;

  G4ThreeVector AdjustSecondaryPosition(G4ThreeVector pos) const;

protected:
  const G4LatticePhysical* theLattice;	// For convenient access by processes

  const G4ParticleDefinition* GetCurrentParticle() const;
  const G4Track* GetCurrentTrack() const { return currentTrack; }
  const G4VPhysicalVolume* GetCurrentVolume() const { return currentVolume; }

  G4int GetCurrentValley() const { return GetValleyIndex(currentTrack); }
  G4int fPhysicsModelID;

private:
  const G4Track* currentTrack;		// For use by Start/EndTracking
  const G4VPhysicalVolume* currentVolume;

  G4AffineTransform fLocalToGlobal;	// For converting pos and momentum
  G4AffineTransform fGlobalToLocal;
};

#endif	/* G4CMPProcessUtils_hh */

/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPBoundaryUtils.hh
/// \brief Definition of the G4CMPBoundaryUtils class
///   Provides useful general functions for use by boundary processes,
///   which do not have a common base class for inheritance.
///
///   Use via multiple inheritance with concrete boundary classes.
//
// $Id$
//
// 20160904  Add electrode pattern handling
// 20160906  Make most functions const, provide casting function for matTable
// 20170525  Add default "rule of five" copy/move operators (no owned pointers)
// 20170627  Make electrode pointer non-const, for initialization
// 20170713  Add registry to keep track of missing-surface warnings
// 20171215  Change 'CheckStepStatus()' to 'IsBoundaryStep()', add function
//	     to validate step trajectory to boundary.
// 20250927  Add overloadable function to kill track when max-reflections.
// 20251204  G4CMP-511 -- Create parallel Lambertian reflection code for charges.
// 20251210  G4CMP-518 -- Make PhononVelocityIsInward() generic.

#ifndef G4CMPBoundaryUtils_hh
#define G4CMPBoundaryUtils_hh 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include <map>
#include <utility>

class G4CMPProcessUtils;
class G4CMPSurfaceProperty;
class G4CMPVElectrodePattern;
class G4LatticePhysical;
class G4MaterialPropertiesTable;
class G4ParticleChange;
class G4Step;
class G4Track;
class G4VPhysicalVolume;
class G4VProcess;


class G4CMPBoundaryUtils {
public:
  G4CMPBoundaryUtils(G4VProcess* process);
  virtual ~G4CMPBoundaryUtils();

  G4CMPBoundaryUtils(const G4CMPBoundaryUtils&) = default;
  G4CMPBoundaryUtils(G4CMPBoundaryUtils&&) = default;
  G4CMPBoundaryUtils& operator=(const G4CMPBoundaryUtils&) = default;
  G4CMPBoundaryUtils& operator=(G4CMPBoundaryUtils&&) = default;
  
  virtual void SetVerboseLevel(G4int vb) { buVerboseLevel = vb; }

  // Check whether this step is at a good boundary for processing
  virtual G4bool IsGoodBoundary(const G4Step& aStep);

  // Check whether end of step is actually on surface of volume
  // "surfacePoint" returns post-step position, or computed surface point
  virtual G4bool CheckStepBoundary(const G4Step& aStep,
				   G4ThreeVector& surfacePoint);

  // Implements PostStepDoIt() in a common way; processes should call through
  virtual void ApplyBoundaryAction(const G4Track& aTrack, const G4Step& aStep,
				   G4ParticleChange& aParticleChange);

  // Decide and apply different surface actions; subclasses may override
  virtual G4bool AbsorbTrack(const G4Track& aTrack, const G4Step& aStep) const;
  virtual void DoAbsorption(const G4Track& aTrack, const G4Step& aStep,
			    G4ParticleChange& aParticleChange);

  virtual G4bool ReflectTrack(const G4Track& aTrack, const G4Step& aStep) const;
  virtual void DoReflection(const G4Track& aTrack, const G4Step& aStep,
			    G4ParticleChange& aParticleChange);

  virtual G4bool MaximumReflections(const G4Track& aTrack) const;
  virtual void DoFinalReflection(const G4Track& aTrack, const G4Step& aStep,
				 G4ParticleChange& aParticleChange);
  // DriftBoundaryProcess should override above to handle recombination

  // Minimal action to simply remove track; no energy transfer, no secondaries.
  virtual void DoSimpleKill(const G4Track& aTrack, const G4Step& aStep,
			    G4ParticleChange& aParticleChange);

  // NOTE:  Transmission is called only if absorption, reflection both fail
  virtual void DoTransmission(const G4Track& aTrack, const G4Step& aStep,
			      G4ParticleChange& aParticleChange);

  // Phonons reflect difusively from surfaces.
  virtual G4ThreeVector LambertianReflection(const G4LatticePhysical* theLattice,
                                    const G4ThreeVector& surfNorm, G4int mode);
  virtual G4ThreeVector LambertianReflection(const G4LatticePhysical* theLattice,
                                    const G4ThreeVector& surfNorm, G4int mode,
                                    const G4ThreeVector& surfPoint);
  virtual G4ThreeVector GetLambertianVector(const G4ThreeVector& surfNorm);

  // Test that a phonon's wave vector relates to an inward velocity.
  // waveVector, surfNorm, and surfacePos need to be in global coordinates
  virtual G4bool VelocityIsInward(const G4LatticePhysical* lattice, G4int mode,
                                const G4ThreeVector& waveVector,
                                const G4ThreeVector& surfNorm);
  virtual G4bool VelocityIsInward(const G4LatticePhysical* lattice, G4int mode,
                                const G4ThreeVector& waveVector,
                                const G4ThreeVector& surfNorm,
                                const G4ThreeVector& surfacePos);

protected:
  G4bool IsBounaryStep(const G4Step& aStep);
  G4bool GetBoundingVolumes(const G4Step& aStep);
  G4bool GetSurfaceProperty(const G4Step& aStep);

  // Does const-casting of matTable for access
  G4double GetMaterialProperty(const G4String& key) const;

private:
  G4int buVerboseLevel;			// For local use; name avoids collisions
  G4String procName;
  G4CMPProcessUtils* procUtils;		// For access to lattice, track info

protected:
  G4double kCarTolerance;		// Allowed nearness to surface
  G4int maximumReflections;		// Limit on track reflections
  G4VPhysicalVolume* prePV;		// Volumes on each side of boundary
  G4VPhysicalVolume* postPV;
  G4CMPSurfaceProperty* surfProp;	// Surface property with G4CMP data
  G4MaterialPropertiesTable* matTable;	// Phonon- or charge-specific parameters
  G4CMPVElectrodePattern* electrode;	// Patterned electrode for absorption

  // Flag whether a given PV pair has a defined surface property or not
  typedef std::pair<G4VPhysicalVolume*,G4VPhysicalVolume*> BoundaryPV;
  std::map<BoundaryPV, G4bool> hasSurface;
};

#endif	/* G4CMPBoundaryUtils_hh */

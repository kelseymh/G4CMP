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

#ifndef G4CMPBoundaryUtils_hh
#define G4CMPBoundaryUtils_hh 1

#include "globals.hh"
#include "G4ThreeVector.hh"

class G4CMPProcessUtils;
class G4CMPSurfaceProperty;
class G4CMPVElectrodePattern;
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

  // Implements PostStepDoIt() in a common way; processes should call through
  virtual G4bool ApplyBoundaryAction(const G4Track& aTrack, const G4Step& aStep,
				     G4ParticleChange& aParticleChange);
  // NOTE:  If step condition tests fail, function will return false, and
  //	    processes should call through to G4VDiscreteProcess:PostStepDoIt()

  // Decide and apply different surface actions; subclasses may override
  virtual G4bool AbsorbTrack(const G4Track& aTrack, const G4Step& aStep);
  virtual void DoAbsorption(const G4Track& aTrack, const G4Step& aStep,
			    G4ParticleChange& aParticleChange);

  virtual G4bool ReflectTrack(const G4Track& aTrack, const G4Step& aStep);
  virtual void DoReflection(const G4Track& aTrack, const G4Step& aStep,
			    G4ParticleChange& aParticleChange);

  virtual G4bool MaximumReflections(const G4Track& aTrack);
  virtual void DoSimpleKill(const G4Track& aTrack, const G4Step& aStep,
			    G4ParticleChange& aParticleChange);

  // NOTE:  Transmission is called only if absorption, reflection both fail
  virtual void DoTransmission(const G4Track& aTrack, const G4Step& aStep,
			      G4ParticleChange& aParticleChange);

protected:
  // Initialize volumes, surface properties, etc.
  G4bool LoadDataForStep(const G4Step& aStep);
  G4bool CheckStepStatus(const G4Step& aStep);
  G4bool GetBoundingVolumes(const G4Step& aStep);
  G4bool GetSurfaceProperty(const G4Step& aStep);

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
  const G4CMPVElectrodePattern* electrode; // Patterned electrode for absorption
};

#endif	/* G4CMPBoundaryUtils_hh */

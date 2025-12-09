/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPDriftBoundaryProcess.hh
/// \brief Subclass of G4CMPVBoundaryProcess to implement interaction of
///        charges (electrons and holes) with crystal surfaces.
//
// $Id$
//
// 20140313  Introduce multiple inheritance from G4CMPProcessUtils
// 20140331  Inherit from G4CMPVDriftProcess to get subtype enforcement
// 20150212  Remove file IO. Use sensitive detectors instead
// 20150304  Change to generic G4CMPVDrifBoundaryProcess and
//           utilize specific G4CMPDrift{Electron,Hole}BoundaryProcess
// 20150420  Replace MFP with GPIL to suppress unnecessary verbosity.
// 20150529  Add DoReflection() function, so electron can overload
// 20150603  Add parameter to count number of reflections by track
// 20160215  Reunify boundary processes at least until we remove CDMS-specific
//           stuff.
// 20160415  Refactor to make PostStepDoIt() very generic. Derived classes
//           should only need to implement absorb/reflect/transmit functions.
// 20160903  Add inheritance from G4CMPBoundaryUtils, remove redundant functions
// 20160906  Follow constness of G4CMPBoundaryUtils
// 20170731  Split electron, hole reflection into utility functions
// 20170802  Add EnergyPartition to handle phonon production
// 20250927  Add override version of new DoFinalReflection(), to support
//           proper recombination.
// 20251013  Add functions for specular and diffuse electron reflection.
// 20251204  G4CMP-511 -- Create parallel Lambertian reflection code for charges.

#ifndef G4CMPDriftBoundaryProcess_h
#define G4CMPDriftBoundaryProcess_h 1

#include "G4CMPVDriftProcess.hh"
#include "G4CMPBoundaryUtils.hh"

class G4CMPEnergyPartition;


class G4CMPDriftBoundaryProcess : public G4CMPVDriftProcess,
				  public G4CMPBoundaryUtils {
public:
  G4CMPDriftBoundaryProcess(const G4String& name = "G4CMPChargeBoundary");
  virtual ~G4CMPDriftBoundaryProcess();

  virtual G4double PostStepGetPhysicalInteractionLength(const G4Track& track,
                                                   G4double previousStepSize,
                                                   G4ForceCondition* condition);

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

protected:
  virtual G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*);

  // Decide and apply different surface actions; subclasses may override
  virtual G4bool AbsorbTrack(const G4Track& aTrack, const G4Step& aStep) const;

  virtual void DoAbsorption(const G4Track& aTrack, const G4Step& aStep,
			    G4ParticleChange& aParticleChange);

  virtual void DoReflection(const G4Track& aTrack,const G4Step& aStep,
			    G4ParticleChange& aParticleChange);

  virtual void DoReflectionElectron(const G4Track& aTrack,const G4Step& aStep,
				    G4ParticleChange& aParticleChange);

  virtual G4ThreeVector DoSpecularElectron(const G4ThreeVector& inDir,
					   const G4ThreeVector& surfNorm,
					   const G4ThreeVector& surfPos) const;

  virtual G4ThreeVector DoDiffuseElectron(const G4ThreeVector& surfNorm,
					  const G4ThreeVector& surfPos) const;

  virtual void DoReflectionHole(const G4Track& aTrack,const G4Step& aStep,
				G4ParticleChange& aParticleChange);

  // Called when maximum bounces have been recorded; does recombination
  virtual void DoFinalReflection(const G4Track& aTrack,const G4Step& aStep,
				 G4ParticleChange& aParticleChange);

  // Phonons reflect difusively from surfaces.
  G4ThreeVector LambertianReflection(const G4LatticePhysical* theLattice,
                                    const G4ThreeVector& surfNorm, G4int valley) override;
  G4ThreeVector LambertianReflection(const G4LatticePhysical* theLattice,
                                    const G4ThreeVector& surfNorm, G4int valley,
                                    const G4ThreeVector& surfPoint) override;

  // Test that a charge's wave vector relates to an inward velocity.
  // waveVector, surfNorm, and surfacePos need to be in global coordinates
  virtual G4bool ChargeVelocityIsInward(const G4LatticePhysical* lattice, G4int valley,
                                const G4ThreeVector& waveVector,
                                const G4ThreeVector& surfNorm);
  virtual G4bool ChargeVelocityIsInward(const G4LatticePhysical* lattice, G4int valley,
                                const G4ThreeVector& waveVector,
                                const G4ThreeVector& surfNorm,
                                const G4ThreeVector& surfacePos);

private:
  G4CMPEnergyPartition* partitioner;

  // No copying/moving
  G4CMPDriftBoundaryProcess(G4CMPDriftBoundaryProcess&);
  G4CMPDriftBoundaryProcess(G4CMPDriftBoundaryProcess&&);
  G4CMPDriftBoundaryProcess& operator=(const G4CMPDriftBoundaryProcess&);
  G4CMPDriftBoundaryProcess& operator=(const G4CMPDriftBoundaryProcess&&);
};

#endif	/* G4CMPDriftBoundaryProcess_h */

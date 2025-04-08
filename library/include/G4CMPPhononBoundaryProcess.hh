/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPPhononBoundaryProcess.hh
/// \brief Definition of the G4CMPPhononBoundaryProcess class
//
// $Id$
//
// 20160903  Add inheritance from G4CMPBoundaryUtils, remove redundant functions
// 20160906  Follow constness of G4CMPBoundaryUtils
// 20181010  J. Singh -- Use new G4CMPAnharmonicDecay for boundary decays
// 20181011  M. Kelsey -- Add LoadDataForTrack() to initialize decay utility.
// 20220906  M. Kelsey -- Encapsulate specular reflection in function.
// 20250204  N. Tenpas -- Support reflection displacement search at hard corners.
// 20250325  N. Tenpas -- Add support for macro commands to set step size and limit.

#ifndef G4CMPPhononBoundaryProcess_h
#define G4CMPPhononBoundaryProcess_h 1

#include "G4VPhononProcess.hh"
#include "G4CMPBoundaryUtils.hh"

class G4CMPAnharmonicDecay;

class G4CMPPhononBoundaryProcess : public G4VPhononProcess,
				   public G4CMPBoundaryUtils {
public:
  G4CMPPhononBoundaryProcess(const G4String& processName="G4CMPPhononBoundary");

  virtual ~G4CMPPhononBoundaryProcess();

  // Configure for current track including AnharmonicDecay utility
  virtual void LoadDataForTrack(const G4Track* track);

  virtual G4double PostStepGetPhysicalInteractionLength(const G4Track& track,
                                                G4double previousStepSize,
                                                G4ForceCondition* condition);

  virtual G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                          const G4Step& aStep);

protected:
  virtual G4double GetMeanFreePath(const G4Track& aTrack,
                                   G4double prevStepLength,
                                   G4ForceCondition* condition);

  // Apply phonon-specific conditions, after calling through to base
  virtual G4bool AbsorbTrack(const G4Track& aTrack, const G4Step& aStep) const;

  // Determine if we are in the special displaced-reflection case
  G4bool BoundaryToBoundaryStep(const G4Step& aStep);

  virtual void DoReflection(const G4Track& aTrack, const G4Step& aStep,
			                      G4ParticleChange& aParticleChange);

  G4ThreeVector GetReflectedVector(const G4ThreeVector& waveVector, 
				                           G4ThreeVector& surfNorm, G4int mode,
                                   G4ThreeVector& surfacePoint);

  G4ThreeVector GetLambertianVector(const G4ThreeVector& surfNorm,
				                            G4int mode) const;


  // Efficiently find direction with min distance to surface
  void OptimizeSurfaceAdjustAngle(const G4VSolid* solid,
                                  const G4ThreeVector& stepLocalPos,
                                  G4double& theta0, G4double& phi0,
                                  const G4int angOption,
                                  const G4double minDist) const;

  // Modifies stepLocalPos in place
  void AdjustToClosestSurfacePoint(const G4VSolid* solid,
                                   G4ThreeVector& stepLocalPos) const;

  void AdjustToEdgePosition(const G4VSolid* solid, const G4ThreeVector& kTan,
                            G4ThreeVector& stepLocalPos) const;

  // Modifies kTan and newNorm in place
  void ReflectAgainstEdge(const G4VSolid* solid, G4ThreeVector& kTan,
                          const G4ThreeVector& stepLocalPos, G4ThreeVector& newNorm) const;

private:
  G4CMPAnharmonicDecay* anharmonicDecay;
  G4double stepSize;
  G4int nStepLimit;

  // hide assignment operator as private
  G4CMPPhononBoundaryProcess(G4CMPPhononBoundaryProcess&);
  G4CMPPhononBoundaryProcess(G4CMPPhononBoundaryProcess&&);
  G4CMPPhononBoundaryProcess& operator=(const G4CMPPhononBoundaryProcess&);
  G4CMPPhononBoundaryProcess& operator=(const G4CMPPhononBoundaryProcess&&);
};

#endif	/* G4CMPPhononBoundaryProcess_h */

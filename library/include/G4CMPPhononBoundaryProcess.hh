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
// 20250422  N. Tenpas -- Add position argument to GetLambertianVector.
// 20250423  G4CMP-468 -- Add wrapper function for updating navigator.
// 20250423  G4CMP-468 -- Move GetLambertianVector to G4CMPUtils.
// 20250424  G4CMP-465 -- Add G4CMPSolidUtils object for custom solid functions.

#ifndef G4CMPPhononBoundaryProcess_h
#define G4CMPPhononBoundaryProcess_h 1

#include "G4VPhononProcess.hh"
#include "G4CMPBoundaryUtils.hh"
#include "G4CMPParticleChangeForPhonon.hh"

class G4CMPAnharmonicDecay;
class G4CMPSolidUtils;

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

  virtual void DoReflection(const G4Track& aTrack, const G4Step& aStep,
			                      G4ParticleChange& aParticleChange);

  G4ThreeVector GetSpecularVector(const G4ThreeVector& waveVector,
                                  G4ThreeVector& surfNorm, G4int mode,
                                  G4ThreeVector& surfacePoint);

  // Update navigator volume when position is changed
  void UpdateNavigatorVolume(const G4Step&, const G4ThreeVector& position,
                             const G4ThreeVector& vDir) const;

private:
  G4CMPAnharmonicDecay* anharmonicDecay;
  G4CMPParticleChangeForPhonon phParticleChange;
  G4double stepSize;
  G4int nStepLimit;

  // hide assignment operator as private
  G4CMPPhononBoundaryProcess(G4CMPPhononBoundaryProcess&);
  G4CMPPhononBoundaryProcess(G4CMPPhononBoundaryProcess&&);
  G4CMPPhononBoundaryProcess& operator=(const G4CMPPhononBoundaryProcess&);
  G4CMPPhononBoundaryProcess& operator=(const G4CMPPhononBoundaryProcess&&);
};

#endif	/* G4CMPPhononBoundaryProcess_h */

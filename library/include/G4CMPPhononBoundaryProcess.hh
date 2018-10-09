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

#ifndef G4CMPPhononBoundaryProcess_h
#define G4CMPPhononBoundaryProcess_h 1

#include "G4VPhononProcess.hh"
#include "G4CMPBoundaryUtils.hh"


class G4CMPPhononBoundaryProcess : public G4VPhononProcess,
				   public G4CMPBoundaryUtils {
public:
  G4CMPPhononBoundaryProcess(const G4String& processName="G4CMPPhononBoundary");

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

private:
  // hide assignment operator as private
  G4CMPPhononBoundaryProcess(G4CMPPhononBoundaryProcess&);
  G4CMPPhononBoundaryProcess(G4CMPPhononBoundaryProcess&&);
  G4CMPPhononBoundaryProcess& operator=(const G4CMPPhononBoundaryProcess&);
  G4CMPPhononBoundaryProcess& operator=(const G4CMPPhononBoundaryProcess&&);
	const G4double BoundaryAnharmonicProb(const G4Track&);
	const G4double BoundarySpecularProb(const G4double);
	const G4double BoundaryLambertianProb(const G4double);
	G4ThreeVector GetLambertianVector(G4ThreeVector);
};

#endif	/* G4CMPPhononBoundaryProcess_h */

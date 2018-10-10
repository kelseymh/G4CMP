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

  // Compute probabilities for different reflection actions
  G4double AnharmonicProb(G4double) const;
  G4double SpecularProb(G4double) const;
  G4double DiffuseProb(G4double) const;
  G4ThreeVector GetLambertianVector(const G4ThreeVector&, G4int) const;

private:
  G4CMPAnharmonicDecay* anharmonicDecay;

  // hide assignment operator as private
  G4CMPPhononBoundaryProcess(G4CMPPhononBoundaryProcess&);
  G4CMPPhononBoundaryProcess(G4CMPPhononBoundaryProcess&&);
  G4CMPPhononBoundaryProcess& operator=(const G4CMPPhononBoundaryProcess&);
  G4CMPPhononBoundaryProcess& operator=(const G4CMPPhononBoundaryProcess&&);
};

#endif	/* G4CMPPhononBoundaryProcess_h */

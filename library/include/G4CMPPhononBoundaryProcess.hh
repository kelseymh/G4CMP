/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPPhononBoundaryProcess.hh
/// \brief Definition of the G4CMPPhononBoundaryProcess class
//
// $Id$
//
#ifndef G4CMPPhononBoundaryProcess_h
#define G4CMPPhononBoundaryProcess_h 1

#include "G4CMPVBoundaryProcessCommon.hh"
#include "G4VPhononProcess.hh"

class G4CMPTrackInformation;
class G4SurfaceProperty;

class G4CMPPhononBoundaryProcess : public G4VPhononProcess {
public:
  G4CMPPhononBoundaryProcess(const G4String& processName =
                             "G4CMPPhononBoundary");
  virtual G4double PostStepGetPhysicalInteractionLength(const G4Track& track,
                                                G4double previousStepSize,
                                                G4ForceCondition* condition);
  
  virtual G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                          const G4Step& aStep);
  // NOTE: This function and any dervied must call back to base class implementation!
  virtual void LoadDataForTrack(const G4Track* track);
  
protected:
  virtual G4double GetMeanFreePath(const G4Track& aTrack,
                                   G4double prevStepLength,
                                   G4ForceCondition* condition);

  // Decide and apply different surface actions; subclasses may override
  virtual G4bool AbsorbTrack(const G4Track& aTrack,
                             const G4Step& aStep,
                             const G4SurfaceProperty* surfProp);

  virtual G4VParticleChange* DoAbsorption(const G4Track& aTrack,
                                          const G4Step& aStep,
                                          const G4SurfaceProperty* surfProp);

  virtual G4bool ReflectTrack(const G4Track& aTrack,
                              const G4Step& aStep,
                              const G4SurfaceProperty* surfProp);

  virtual G4VParticleChange* DoReflection(const G4Track& aTrack,
                                          const G4Step& aStep,
                                          const G4SurfaceProperty* surfProp);

  virtual G4VParticleChange* DoTransmission(const G4Track& aTrack,
                                            const G4Step& aStep,
                                            const G4SurfaceProperty* surfProp);

  G4bool ReflectionIsGood(G4int polarization, G4ThreeVector waveVector,
                          G4ThreeVector surfNorm);

  G4double kCarTolerance;
  G4double maxNumReflections;

private:
  // hide assignment operator as private
  G4CMPPhononBoundaryProcess(G4CMPPhononBoundaryProcess&);
  G4CMPPhononBoundaryProcess(G4CMPPhononBoundaryProcess&&);
  G4CMPPhononBoundaryProcess& operator=(const G4CMPPhononBoundaryProcess&);
  G4CMPPhononBoundaryProcess& operator=(const G4CMPPhononBoundaryProcess&&);
};

#endif











/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4PhononReflection.hh
/// \brief Definition of the G4PhononReflection class
//
// $Id$
//
#ifndef G4PhononReflection_h
#define G4PhononReflection_h 1

#include "G4VPhononProcess.hh"

class G4CMPTrackInformation;
class G4CMPSurfaceProperty;

class G4PhononReflection : public G4VPhononProcess {
public:
  G4PhononReflection(const G4String& processName ="phononReflection" );
  
  virtual G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                          const G4Step& aStep);
  
protected:
  virtual G4double GetMeanFreePath(const G4Track& aTrack,
                                   G4double prevStepLength,
                                   G4ForceCondition* condition);
  // Decide and apply different surface actions; subclasses may override
  virtual G4bool AbsorbTrack(const G4Step& aStep);
  virtual G4VParticleChange* DoAbsorption(const G4Step& aStep);

  virtual G4bool ReflectTrack(const G4Step& aStep);
  virtual G4VParticleChange* DoReflection(const G4Step& aStep);

  G4bool ReflectionIsGood(G4int polarization);

  G4double kCarTolerance;
  G4double absProb;
  G4double reflProb;
  G4double specProb;
  G4double absMinK;
  G4double maxRefl;
  G4CMPTrackInformation* trackInfo;
  G4ThreeVector waveVector;
  G4ThreeVector surfNorm;

private:
  // hide assignment operator as private
  G4PhononReflection(G4PhononReflection&);
  G4PhononReflection& operator=(const G4PhononReflection& right);
};

#endif











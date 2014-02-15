#ifndef G4CMPTimeStepper_h
#define G4CMPTimeStepper_h 1

#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "G4AffineTransform.hh"

class G4CMPTimeStepper : public G4VDiscreteProcess {
public:
  G4CMPTimeStepper(const G4String& processName="G4CMPTimeStepper");

  virtual ~G4CMPTimeStepper();

  virtual G4bool IsApplicable(const G4ParticleDefinition &pd);

  virtual G4double 
  PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
				       G4double prevStepSize,
				       G4ForceCondition* condition);

  virtual G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
					  const G4Step& aStep);

protected:
  // Compute dt_e, dt_h and valley rotations at current location
  void ComputeTimeSteps(const G4Track& aTrack);

  virtual G4double GetMeanFreePath(const G4Track& aTrack,
				   G4double prevStepSize,
				   G4ForceCondition* condition) {
    return PostStepGetPhysicalInteractionLength(aTrack, prevStepSize,
						condition);
  }

private:
  //hide assignment operator
  G4CMPTimeStepper(G4CMPTimeStepper&);
  G4CMPTimeStepper& operator=(const G4CMPTimeStepper& right);

  G4AffineTransform normalToValley;
  G4AffineTransform valleyToNormal;
  G4double dt_e; //Time step for electrons
  G4double dt_h; //Time step for holes

  const G4double velLong;
  const G4double me;
  const G4double mc_e;
  const G4double l0_e;
  const G4double ksound_e;
  const G4double mc_h;
  const G4double l0_h;
  const G4double ksound_h;
  const G4ThreeVector T;
};

#endif

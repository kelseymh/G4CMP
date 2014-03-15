// $Id$
//
// 20140313  Introduce multiple inheritance from G4CMPProcessUtils
//	     Add wrapper function to compute individual time steps

#ifndef G4CMPTimeStepper_h
#define G4CMPTimeStepper_h 1

#include "globals.hh"
#include "G4CMPVDriftProcess.hh"
#include "G4AffineTransform.hh"


class G4CMPTimeStepper : public G4CMPVDriftProcess {
public:
  G4CMPTimeStepper();
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
  G4double TimeStepInField(G4double Efield, G4double coeff, G4double l0) const;

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
};

#endif

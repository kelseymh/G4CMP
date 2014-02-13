#ifndef G4CMPTimeStepper_h
#define G4CMPTimeStepper_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VProcess.hh"
#include "G4AffineTransform.hh"

class G4CMPTimeStepper : public G4VProcess
{
public:
  G4CMPTimeStepper(const G4String& processName="G4CMPTimeStepper");

  virtual ~G4CMPTimeStepper();

  virtual G4double PostStepGetPhysicalInteractionLength(
							const G4Track& aTrack,
							G4double prevStepSize,
							G4ForceCondition* condition
							);

  virtual G4VParticleChange* PostStepDoIt(
					  const G4Track& aTrack,
					  const G4Step& aStep
					  );

  //no-op in atRest GPIL
  virtual G4double AtRestGetPhysicalInteractionLength(
						      const G4Track&,
						      G4ForceCondition*
                  ){ return -1.0;}

  //no-op in AtRestDoIt
  virtual G4VParticleChange* AtRestDoIt(
					const G4Track&,
					const G4Step&
          ){ return 0;}

  //no-op in alongStepGPIL
  virtual G4double AlongStepGetPhysicalInteractionLength(
						      const G4Track&,
						      G4double,
						      G4double,
						      G4double&,
						      G4GPILSelection*
                  ){ return -1.0;}

  //no-op in alongStepDoIt
  virtual G4VParticleChange* AlongStepDoIt(
					const G4Track&,
					const G4Step&
          ){ return 0;}

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

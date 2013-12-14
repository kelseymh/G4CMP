#ifndef ObliqueEFieldProcess_h
#define ObliqueEFieldProcess_h 1

#include "G4VDiscreteProcess.hh"
#include "G4AffineTransform.hh"

class ObliqueEFieldProcess : public G4VDiscreteProcess
{

public:

  G4AffineTransform normalToValley;
  G4AffineTransform valleyToNormal;

  ObliqueEFieldProcess(const G4String& processName = "EfieldProcess");

  virtual ~ObliqueEFieldProcess();

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

  virtual bool IsApplicable(const G4ParticleDefinition&);

protected:

  virtual G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*);

private:

  ObliqueEFieldProcess(ObliqueEFieldProcess&);
  ObliqueEFieldProcess& operator=(const ObliqueEFieldProcess& right);

};
 #endif

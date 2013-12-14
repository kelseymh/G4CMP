#ifndef EFieldProcess_h
#define EFieldProcess_h 1

#include "G4VDiscreteProcess.hh"

class EFieldProcess : public G4VDiscreteProcess
{

public:


  EFieldProcess(const G4String& processName = "EfieldProcess");

  virtual ~EFieldProcess();

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

  virtual bool IsApplicable(const G4ParticleDefinition&);

protected:

  virtual G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*);

private:

  EFieldProcess(EFieldProcess&);
  EFieldProcess& operator=(const EFieldProcess& right);

};
 #endif

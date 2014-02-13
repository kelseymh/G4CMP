#ifndef G4CMPDriftBoundaryProcess_h
#define G4CMPDriftBoundaryProcess_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VDiscreteProcess.hh"

class G4CMPDriftBoundaryProcess : public G4VDiscreteProcess
{
public:

  G4CMPDriftBoundaryProcess(const G4String& processName = "G4CMPDriftBoundaryProcess");

  virtual ~G4CMPDriftBoundaryProcess();

  virtual G4VParticleChange* PostStepDoIt(
                const G4Track&, const G4Step& );

  virtual G4bool IsApplicable(const G4ParticleDefinition&);

protected:
  virtual G4double GetMeanFreePath(
                const G4Track&, G4double, G4ForceCondition* );

private:
  G4CMPDriftBoundaryProcess(G4CMPDriftBoundaryProcess&);
  G4CMPDriftBoundaryProcess& operator = (const G4CMPDriftBoundaryProcess& right);

  G4double kCarTolerance;

};

#endif

#ifndef DriftingCarrierBoundaryProcess_h
#define DriftingCarrierBoundaryProcess_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VDiscreteProcess.hh"

class DriftingCarrierBoundaryProcess : public G4VDiscreteProcess
{
public:

  DriftingCarrierBoundaryProcess(const G4String& processName = "DriftingCarrierBoundaryProcess");

  virtual ~DriftingCarrierBoundaryProcess();

  virtual G4VParticleChange* PostStepDoIt(
                const G4Track&, const G4Step& );

  virtual G4bool IsApplicable(const G4ParticleDefinition&);

protected:
  virtual G4double GetMeanFreePath(
                const G4Track&, G4double, G4ForceCondition* );

private:
  DriftingCarrierBoundaryProcess(DriftingCarrierBoundaryProcess&);
  DriftingCarrierBoundaryProcess& operator = (const DriftingCarrierBoundaryProcess& right);

  G4double kCarTolerance;

};

#endif

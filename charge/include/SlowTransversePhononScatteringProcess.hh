#ifndef SlowTransversePhononScatteringProcess_h
#define SlowTransversePhononScatteringProcess_h 1
#define ST 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "LatticeManager2.hh"

class SlowTransversePhononScatteringProcess : public G4VDiscreteProcess 
{
  public:

     SlowTransversePhononScatteringProcess(const G4String& processName ="SlowTransversePhononScattering" );

  //     SlowTransversePhononScatteringProcess();

     virtual ~SlowTransversePhononScatteringProcess();

     virtual G4VParticleChange* PostStepDoIt(
                const G4Track&, const G4Step& );
 
     virtual G4bool IsApplicable(const G4ParticleDefinition&);
                           
  protected:

     virtual G4double GetMeanFreePath(
                const G4Track&, G4double, G4ForceCondition* );

  private:
  LatticeManager2 LM;
  
  // hide assignment operator as private 
     SlowTransversePhononScatteringProcess(SlowTransversePhononScatteringProcess&);
     SlowTransversePhononScatteringProcess& operator=(const SlowTransversePhononScatteringProcess& right);

};

#endif











#ifndef TransversePhononScatteringProcess_h
#define TransversePhononScatteringProcess_h 1
#define TRANS 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "LatticeManager2.hh"

class TransversePhononScatteringProcess : public G4VDiscreteProcess 
{
  public:

     TransversePhononScatteringProcess(const G4String& processName ="TransversePhononScattering" );

  //     TransversePhononScatteringProcess();

     virtual ~TransversePhononScatteringProcess();

     virtual G4VParticleChange* PostStepDoIt(
                const G4Track&, const G4Step& );
 
     virtual G4bool IsApplicable(const G4ParticleDefinition&);
                           
  protected:

     virtual G4double GetMeanFreePath(
                const G4Track&, G4double, G4ForceCondition* );

  private:
  LatticeManager2 LM;
  
  // hide assignment operator as private 
     TransversePhononScatteringProcess(TransversePhononScatteringProcess&);
     TransversePhononScatteringProcess& operator=(const TransversePhononScatteringProcess& right);

};

#endif











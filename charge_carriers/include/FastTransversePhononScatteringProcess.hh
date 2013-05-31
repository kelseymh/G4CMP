#ifndef FastTransversePhononScatteringProcess_h
#define FastTransversePhononScatteringProcess_h 1
#define FT 2

#include "G4ios.hh"
#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "LatticeManager2.hh"

class FastTransversePhononScatteringProcess : public G4VDiscreteProcess 
{
  public:

     FastTransversePhononScatteringProcess(const G4String& processName ="FastTransversePhononScattering" );

  //     FastTransversePhononScatteringProcess();

     virtual ~FastTransversePhononScatteringProcess();

     virtual G4VParticleChange* PostStepDoIt(
                const G4Track&, const G4Step& );
 
     virtual G4bool IsApplicable(const G4ParticleDefinition&);
                           
  protected:

     virtual G4double GetMeanFreePath(
                const G4Track&, G4double, G4ForceCondition* );

  private:
  LatticeManager2 LM;
  
  // hide assignment operator as private 
     FastTransversePhononScatteringProcess(FastTransversePhononScatteringProcess&);
     FastTransversePhononScatteringProcess& operator=(const FastTransversePhononScatteringProcess& right);

};

#endif











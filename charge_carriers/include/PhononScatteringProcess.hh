#ifndef PhononScatteringProcess_h
#define PhononScatteringProcess_h 1
#define LON 0

#include "G4ios.hh"
#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "LatticeManager2.hh"

class PhononScatteringProcess : public G4VDiscreteProcess 
{
  public:


     PhononScatteringProcess(const G4String& processName ="PhononScattering" );

     virtual ~PhononScatteringProcess();

     virtual G4VParticleChange* PostStepDoIt(
                const G4Track&, const G4Step& );
 
     virtual G4bool IsApplicable(const G4ParticleDefinition&);
                           
  protected:

     virtual G4double GetMeanFreePath(
                const G4Track&, G4double, G4ForceCondition* );

  private:
  
  // hide assignment operator as private 
     PhononScatteringProcess(PhononScatteringProcess&);
     PhononScatteringProcess& operator=(const PhononScatteringProcess& right);


};

#endif











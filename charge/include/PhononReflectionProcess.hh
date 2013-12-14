#ifndef PhononReflectionProcess_h
#define PhononReflectionProcess_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VDiscreteProcess.hh"

class G4Material;

class PhononReflectionProcess : public G4VDiscreteProcess 
{
  public:

     PhononReflectionProcess(const G4String& processName ="PhononReflectionProcess" );

  virtual ~PhononReflectionProcess();
  
  virtual G4VParticleChange* PostStepDoIt(
                const G4Track&, const G4Step& );
 
  virtual G4bool IsApplicable(const G4ParticleDefinition&);

  virtual void BuildPhysicsTable(const G4ParticleDefinition&);
                           
  protected:

  virtual G4double GetMeanFreePath(
                const G4Track&, G4double, G4ForceCondition* );



  private:
  
  // hide assignment operator as private 
     PhononReflectionProcess(PhononReflectionProcess&);
     PhononReflectionProcess& operator=(const PhononReflectionProcess& right);

  private:
 
     G4Material* alminum;
     G4double kCarTolerance;
};

#endif











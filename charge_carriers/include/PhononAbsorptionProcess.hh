#ifndef PhononAbsorptionProcess_h
#define PhononAbsorptionProcess_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VDiscreteProcess.hh"

#include "PhysicalLattice.hh"

class PhononAbsorptionProcess : public G4VDiscreteProcess 
{
  public:


     PhononAbsorptionProcess(const G4String& processName ="PhononAbsorptionProcess" );

     virtual ~PhononAbsorptionProcess();

     virtual G4VParticleChange* PostStepDoIt(
                const G4Track&, const G4Step& );
 
     virtual G4bool IsApplicable(const G4ParticleDefinition&);
                           
  protected:

     virtual G4double GetMeanFreePath(
                const G4Track&, G4double, G4ForceCondition* );



 
  private:
    double BETA, GAMMA, LAMBDA, MU;

  // hide assignment operator as private 
     PhononAbsorptionProcess(PhononAbsorptionProcess&);
     PhononAbsorptionProcess& operator=(const PhononAbsorptionProcess& right);

  // relative probability that anharmonic decay occurs L->L'+T'
  inline double getLTDecayProb(G4double, G4double);
  inline double getTTDecayProb(G4double, G4double);
  inline double makeLDeviation(G4double, G4double);
  inline double makeTTDeviation(G4double, G4double);
  inline double makeTDeviation(G4double, G4double);
  void makeTTSecondaries(const G4Track&);
  void makeLTSecondaries(const G4Track&);



  PhysicalLattice* Lattice;


};

#endif











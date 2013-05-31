
#include "Tst1PhysicsList.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4ios.hh"              

#include "TPhononSlow.hh"
#include "TPhononFast.hh"
#include "LPhonon.hh"

#include "PhononScatteringProcess.hh"
#include "PhononAbsorptionProcess.hh"
#include "PhononReflectionProcess.hh"

#include "DriftingCarrierBoundaryProcess.hh"
#include "LukeScatteringProcess.hh"
#include "ElectronLukeScatteringProcess.hh"
#include "EFieldProcess.hh"
#include "ObliqueEFieldProcess.hh"
#include "TimeStepper.hh"


#include "Tst1EMField.hh"

#include "DriftingElectron.hh"
#include "DriftingHole.hh"
#include "G4StepLimiter.hh"
#include "G4UserSpecialCuts.hh"

#include "G4UserLimits.hh"



Tst1PhysicsList::Tst1PhysicsList():  G4VUserPhysicsList()
{
  G4cout<<"\n\nTst1PhysicsList::constructor: running"<<endl;
  defaultCutValue = DBL_MIN;//100*mm;
  SetVerboseLevel(1);
}

Tst1PhysicsList::~Tst1PhysicsList()
{}

void Tst1PhysicsList::ConstructParticle()
{

  DriftingElectron::Definition();
  DriftingHole::Definition();

  LPhonon::PhononDefinition();
  TPhononFast::PhononDefinition();
  TPhononSlow::PhononDefinition();


}


void Tst1PhysicsList::ConstructProcess()
{

  AddTransportation();  
  /*
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();


    if (particleName == "phonon") {
      pmanager->AddDiscreteProcess(new PhononScatteringProcess());
      pmanager->AddDiscreteProcess(new PhononAbsorptionProcess());
      pmanager->AddDiscreteProcess(new PhononReflectionProcess());
    } 
     
    
    if (particleName == "TPhononSlow") {
      //      G4cout<<"Registering slow transverse scattering..."<<G4endl;
      pmanager->AddDiscreteProcess(new PhononScatteringProcess());      
      pmanager->AddDiscreteProcess(new PhononAbsorptionProcess());
      pmanager->AddDiscreteProcess(new PhononReflectionProcess());
    }
    
    if (particleName == "TPhononFast") {
      //      G4cout<<"Registering fast transverse scattering..."<<G4endl;
      pmanager->AddDiscreteProcess(new PhononScatteringProcess());      
      pmanager->AddDiscreteProcess(new PhononAbsorptionProcess());
      pmanager->AddDiscreteProcess(new PhononReflectionProcess());
    }
    if (particleName == "LPhonon") {
      G4cout<<"Registering longitudinal scattering..."<<G4endl;
      pmanager->AddDiscreteProcess(new PhononScatteringProcess());
      pmanager->AddDiscreteProcess(new PhononAbsorptionProcess());
      pmanager->AddDiscreteProcess(new PhononReflectionProcess());
    } 
  }
  */

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();


    if(particleName=="DriftingHole")
    {
	    TimeStepper* stepperHole = new TimeStepper();
	    pmanager->AddDiscreteProcess(stepperHole);
	    //pmanager->AddDiscreteProcess(new LukeScatteringProcess(stepperHole));
	    pmanager->AddDiscreteProcess(new LukeScatteringProcess());
        //pmanager->AddDiscreteProcess(new EFieldProcess());
	    pmanager->AddDiscreteProcess(new DriftingCarrierBoundaryProcess());
    }   
    if(particleName=="DriftingElectron")
    {
	    TimeStepper* stepperElectron = new TimeStepper();
	    pmanager->AddDiscreteProcess(stepperElectron);
	    pmanager->AddDiscreteProcess(new ElectronLukeScatteringProcess());
	    //pmanager->AddDiscreteProcess(new ObliqueEFieldProcess());
	    pmanager->AddDiscreteProcess(new DriftingCarrierBoundaryProcess());
    }
  }
}

void Tst1PhysicsList::SetCuts()
{
  // These values are used as the default production thresholds
  // for the world volume.
  SetCutsWithDefault();
}



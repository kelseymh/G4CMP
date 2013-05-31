
//Scattering mean free path and amplitude - mean free path is empirical. 
//Scattering amplitude only depends on polarization vector and not on Phonon
//type (i.e. longitudinal, transverse...)

//#include <iostream>
//#include <fstream>
//using std::ifstream;

#include "TransversePhononScatteringProcess.hh"

#include "G4Step.hh"
#include "G4VParticleChange.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "LPhonon.hh"
#include "TPhononFast.hh"
#include "TPhononSlow.hh"


#include "G4TransportationManager.hh"
#include "G4Navigator.hh"

TransversePhononScatteringProcess::TransversePhononScatteringProcess(const G4String& aName)
: G4VDiscreteProcess(aName)
{
  if (verboseLevel>1) {
     G4cout << GetProcessName() << " is created "<< G4endl;
   }
}

TransversePhononScatteringProcess::~TransversePhononScatteringProcess()
{;}

TransversePhononScatteringProcess::TransversePhononScatteringProcess(TransversePhononScatteringProcess& right)
: G4VDiscreteProcess(right)
{;}
 
G4double 
TransversePhononScatteringProcess::GetMeanFreePath( 
       const G4Track& aTrack, G4double /*previousStepSize*/, G4ForceCondition* condition  )
{

  
///////////////////////////////////////////////////////////////////////////
// Physics here !!!
///////////////////////////////////////////////////////////////////////////
  G4double B=3.67e-41*s*s*s;
  G4double h=6.626e-34*m2*kg/s;
  G4double E= aTrack.GetKineticEnergy();
  //G4cout<<"\nTransversePhononScatteringCrossSection::GetMeanFreePath: B E^4/h^4 calculated as: "<<((E/h)*(E/h)*(E/h)*(E/h)*B)<<" Velocity found: "<<aTrack.GetVelocity();
  
  G4double mfp = 1/((E/h)*(E/h)*(E/h)*(E/h)*B)*aTrack.GetVelocity();
 
  //mfp = 1*mm;
   *condition = NotForced;
   return mfp;
}

G4VParticleChange*
  TransversePhononScatteringProcess::PostStepDoIt( const G4Track& aTrack,
                                 const G4Step&  )
{
 

  G4double modeMixer = G4UniformRand();
  G4ThreeVector newDir = G4RandomDirection();
  aParticleChange.Initialize(aTrack);
  G4cout<<"Entering transverse scattering...";

  if(!LatticeManager2::hasLattice(aTrack.GetVolume())){
    aParticleChange.ProposeEnergy(0.);
    aParticleChange.ProposeNonIonizingEnergyDeposit(aTrack.GetKineticEnergy());
    aParticleChange.ProposeTrackStatus(fStopAndKill); 
  } else {
    G4Navigator* theNavigator =
      G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
    newDir = theNavigator->GetGlobalToLocalTransform().TransformAxis(newDir);
    newDir=theNavigator->GetLocalToGlobalTransform().TransformAxis(newDir);
    aParticleChange.ProposeEnergy(0.);
    aParticleChange.ProposeTrackStatus(fStopAndKill);
    aParticleChange.SetNumberOfSecondaries(1);
    if(modeMixer<0.55){
      newDir=LatticeManager2::mapKtoVDir(aTrack.GetVolume(), 1, newDir);
      aParticleChange.AddSecondary(new G4DynamicParticle(TPhononSlow::PhononDefinition(),newDir, aTrack.GetKineticEnergy()));
       G4cout<<"TransverseScatteringProcess::\nScattering to TPhononSlow...";
    }else if(modeMixer<0.6){
      newDir=LatticeManager2::mapKtoVDir(aTrack.GetVolume(), 2, newDir);
       G4cout<<"TransverseScatteringProcess::\nScattering to TPhononFast...";
      aParticleChange.AddSecondary(new G4DynamicParticle(TPhononFast::PhononDefinition(),newDir, aTrack.GetKineticEnergy()));
    } else {
      newDir=LatticeManager2::mapKtoVDir(aTrack.GetVolume(), 0, newDir);
      aParticleChange.AddSecondary(new G4DynamicParticle(LPhonon::PhononDefinition(),newDir, aTrack.GetKineticEnergy()));
       G4cout<<"TransverseScatteringProcess::\nScattering to LPhonon...";
    }
    
  }
  return &aParticleChange;
}



G4bool TransversePhononScatteringProcess::IsApplicable(const G4ParticleDefinition& aPD)
{
   return (true);
}

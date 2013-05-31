
//Scattering mean free path and amplitude - mean free path is empirical. 
//Scattering amplitude only depends on polarization vector and not on Phonon
//type (i.e. longitudinal, transverse...)

//#include <iostream>
//#include <fstream>
//using std::ifstream;

#include "SlowTransversePhononScatteringProcess.hh"

#include "G4Step.hh"
#include "G4VParticleChange.hh"
#include "Randomize.hh"
#include "TPhononFast.hh"
#include "TPhononSlow.hh"


#include "G4TransportationManager.hh"
#include "G4Navigator.hh"

SlowTransversePhononScatteringProcess::SlowTransversePhononScatteringProcess(const G4String& aName)
: G4VDiscreteProcess(aName)
{
  if (verboseLevel>1) {
     G4cout << GetProcessName() << " is created "<< G4endl;
   }
}

SlowTransversePhononScatteringProcess::~SlowTransversePhononScatteringProcess()
{;}

SlowTransversePhononScatteringProcess::SlowTransversePhononScatteringProcess(SlowTransversePhononScatteringProcess& right)
: G4VDiscreteProcess(right)
{;}
 
G4double 
SlowTransversePhononScatteringProcess::GetMeanFreePath( 
       const G4Track& aTrack, G4double /*previousStepSize*/, G4ForceCondition* condition  )
{

  
///////////////////////////////////////////////////////////////////////////
// Physics here !!!
///////////////////////////////////////////////////////////////////////////
  G4double B=3.67e-41*s*s*s;
  G4double h=6.626e-34*m2*kg/s;
  G4double E= aTrack.GetKineticEnergy();
  //G4cout<<"\nSlowTransversePhononScatteringCrossSection::GetMeanFreePath: B E^4/h^4 calculated as: "<<((E/h)*(E/h)*(E/h)*(E/h)*B)<<" Velocity found: "<<aTrack.GetVelocity();
  
  G4double mfp = 1/((E/h)*(E/h)*(E/h)*(E/h)*B)*aTrack.GetVelocity();
  //G4cout<<"\nLongitudinalPhononScatteringProcess::GetMeanFreePath: found mfp "<<mfp;
  //G4double mfp = 1*mm;
   *condition = NotForced;
   return mfp;
}

G4VParticleChange*
  SlowTransversePhononScatteringProcess::PostStepDoIt( const G4Track& aTrack,
                                 const G4Step&  )
{
  aParticleChange.Initialize(aTrack);
  G4ThreeVector currentDir = aTrack.GetMomentumDirection();
  G4Navigator* theNavigator =
     G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  
  G4double phi = pi*G4UniformRand()*rad;
  G4double theta = pi*G4UniformRand()*rad;
  
  G4ThreeVector newDir = currentDir.rotate(G4ThreeVector(0,0,1), theta).rotate(G4ThreeVector(0,1,0), phi);
  newDir = theNavigator->GetGlobalToLocalTransform().TransformAxis(newDir);
  newDir=LatticeManager2::mapKtoVDir(aTrack.GetVolume(), ST, newDir);
  newDir=theNavigator->GetLocalToGlobalTransform().TransformAxis(newDir);
  aParticleChange.ProposeMomentumDirection(newDir);

  return &aParticleChange;
}



G4bool SlowTransversePhononScatteringProcess::IsApplicable(const G4ParticleDefinition& aPD)
{
   return ((&aPD==TPhononFast::PhononDefinition())|| (&aPD==TPhononSlow::PhononDefinition()));
}

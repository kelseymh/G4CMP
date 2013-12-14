#include "EFieldProcess.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"

EFieldProcess::EFieldProcess(const G4String& aName) : G4VDiscreteProcess(aName)
{
  if(verboseLevel>1){
    G4cout<<GetProcessName()<<" is created "<<G4endl;
  }
}

EFieldProcess::~EFieldProcess()
{ ; }

EFieldProcess::EFieldProcess(EFieldProcess& right)
  : G4VDiscreteProcess(right)
{ ; }

G4double EFieldProcess::GetMeanFreePath(const G4Track& aTrack, G4double, G4ForceCondition* condition)
{
  *condition=Forced;
  return DBL_MAX;
}

G4VParticleChange* EFieldProcess::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
  G4double e = 1.0*eplus;
  G4double me=0.35*9.1e-31*kg;
  G4double EField= 20.0*volt/m;
  G4ThreeVector EDir = G4ThreeVector(0.0, 0.0, 1.0).unit();

  //calculate momentum change
  double dv=e*EField*aStep.GetDeltaTime()/me;
  //G4cout << aStep.GetDeltaTime() << G4endl;
  G4ThreeVector deltaV = EDir;
  deltaV.setMag(dv);

  //calculate new momentum
  G4ThreeVector oldV=aTrack.GetMomentumDirection();
  oldV.setMag(aTrack.GetVelocity());
  G4ThreeVector newV;
  newV.set(oldV.getX() + deltaV.getX(),oldV.getY() + deltaV.getY(),oldV.getZ() + deltaV.getZ());
   
  
  G4double dE=0.5*0.35*9.1e-31*((newV.mag()*newV.mag()) - (oldV.mag()*oldV.mag()))/(m/s)/(m/s)/1.602e-19*eV;
  G4double newE=aTrack.GetKineticEnergy()+dE;
  //G4double newE=aTrack.GetKineticEnergy()+0.5*0.35*9.1e-31*((newV.mag()*newV.mag()) - (oldV.mag()*oldV.mag()))/(m/s)/(m/s)/1.602e-19*eV;

  aParticleChange.Initialize(aTrack);

  aParticleChange.ProposeMomentumDirection(newV.unit());
  aParticleChange.ProposeEnergy(newE);

  return &aParticleChange;
}

G4bool EFieldProcess::IsApplicable(const G4ParticleDefinition& aPD){
  return true;
}

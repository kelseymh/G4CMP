#include "TimeStepper_Rob.hh"

//#include "G4TransportationProcessType.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4UserLimits.hh"
#include "G4VParticleChange.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TimeStepper_Rob::TimeStepper_Rob(const G4String& aName)
  : G4VProcess(aName, fGeneral)
{
  //SetProcessSubType(static_cast<int>(USER_SPECIAL_CUTS));

  if(verboseLevel>0)
    {
      G4cout<<GetProcessName()<<" is created "<<G4endl;
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TimeStepper_Rob::~TimeStepper_Rob()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TimeStepper_Rob::TimeStepper_Rob(TimeStepper_Rob& right)
  : G4VProcess(right)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double TimeStepper_Rob::AlongStepGetPhysicalInteractionLength(
							   const G4Track& aTrack,
							   G4double prevStepSize,
                               G4double minStep,
                               G4double& safety,
                               G4GPILSelection* selection)
{

  G4double ProposedStep = DBL_MAX;
  
  /*
  G4UserLimits* pUserLimits = (
			       aTrack.GetVolume()
			       ->GetLogicalVolume()
			       ->GetUserLimits()
			       );
  
  if(pUserLimits)
    {
      G4double tlimit = pUserLimits->GetUserMaxTime(aTrack);
      if(tlimit<DBL_MAX){
	G4double velocity = aTrack.GetVelocity();
	G4double distn = velocity * tlimit;
	ProposedStep=distn;
      }
    }
  */
  G4double tlimit = 2.7571e-2*ns;
  //G4double tlimit = 2e-3*ns;
  
  G4double velocity = aTrack.GetVelocity();
  G4double distn = velocity * tlimit;
  
  ProposedStep=distn;

  return ProposedStep;
}

G4VParticleChange* TimeStepper_Rob::AlongStepDoIt(
					  const G4Track& aTrack,
					  const G4Step& aStep
					  )
{
  aParticleChange.Initialize(aTrack);
  /*
  G4cout<<"\nTimeStepper::PostStepDoIt: Delta time: " 
	<<(aStep.GetDeltaTime() - 0.001*ns)/ns	
	<< " ns"<<G4endl;
  */
  return &aParticleChange;
}

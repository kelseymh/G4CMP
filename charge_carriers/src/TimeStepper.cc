#include "TimeStepper.hh"

//#include "G4TransportationProcessType.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4Field.hh"
#include "Tst1EMField.hh"
#include "G4UserLimits.hh"
#include "G4VParticleChange.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
 #include <fstream>
 #include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TimeStepper::TimeStepper(const G4String& aName)
  : G4VProcess(aName, fGeneral)
{
  //SetProcessSubType(static_cast<int>(USER_SPECIAL_CUTS));

  if(verboseLevel>0)
    {
      G4cout<<GetProcessName()<<" is created "<<G4endl;
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TimeStepper::~TimeStepper()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TimeStepper::TimeStepper(TimeStepper& right)
  : G4VProcess(right)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double TimeStepper::PostStepGetPhysicalInteractionLength(
							   const G4Track& aTrack,
							   G4double prevStepSize,
							   G4ForceCondition* cond)
{
    //set condition to "Forced"
    *cond = NotForced;
  
    G4double velLong=5400*m/s;
    G4double l0Hole = 108e-6*m;
    G4double massFreeElectron=9.1e-31*kg;
    G4double massHole=0.35*massFreeElectron;
    G4double ksound_Hole=1.6306e7/m;
    G4double hbar= 6.5821e-16*eV*s;
    
    G4StepPoint* stepPoint = aTrack.GetStep()->GetPostStepPoint();
    G4FieldManager* fieldMan = G4TransportationManager::GetTransportationManager()->GetFieldManager();
    const G4ElectricField* field = (G4ElectricField*)fieldMan->GetDetectorField();
    
    G4double fieldVal[6];
    G4ThreeVector position3Vec = stepPoint->GetPosition();
    G4double position[4] = {position3Vec[0], position3Vec[1], position3Vec[2], 0};
    field->GetFieldValue(position,  fieldVal);
    //G4cout <<  "Field Value:" << fieldVal[3] <<  " " <<  fieldVal[4] <<  " " <<  fieldVal[5] <<  G4endl;
    G4ThreeVector fieldVec = G4ThreeVector(fieldVal[3],fieldVal[4], fieldVal[5]);
    
    G4double speed = stepPoint->GetVelocity();
    G4double k = speed*massHole/hbar;
    
    G4double dt = (hbar * (ksound_Hole-k))/(fieldVec.mag())/stepPoint->GetCharge();
    //if (dt<0) dt=-dt;
    
    G4double ProposedStep = speed*dt;
    if (ProposedStep < 1e-4*mm) ProposedStep = 1e-4*mm;

    //G4cout <<  "TimeStepper: " <<  ProposedStep <<  G4endl;
    //return ProposedStep;
    return 1e-4*mm;
}

G4VParticleChange* TimeStepper::PostStepDoIt(
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

  /*std::ofstream epositions;
  epositions.open("e-positions.txt", std::ofstream::app);

  epositions << newPosn.getX() << " " << newPosn.getY() << " " << newPosn.getZ() << "\n";
  epositions.close();
  */
  return &aParticleChange;
}

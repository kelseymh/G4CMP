#include "G4CMPTimeStepper.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4Field.hh"
#include "G4UserLimits.hh"
#include "G4VParticleChange.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4CMPDriftHole.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPValleyTrackMap.hh"
#include "G4VDiscreteProcess.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include <fstream>
#include <iostream>
#include <math.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CMPTimeStepper::G4CMPTimeStepper(const G4String& aName)
  : G4VProcess(aName, fGeneral)
{
  //SetProcessSubType(static_cast<int>(USER_SPECIAL_CUTS));

  if(verboseLevel>0)
    {
      G4cout<<GetProcessName()<<" is created "<<G4endl;
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4CMPTimeStepper::~G4CMPTimeStepper()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CMPTimeStepper::G4CMPTimeStepper(G4CMPTimeStepper& right)
  : G4VProcess(right)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4CMPTimeStepper::PostStepGetPhysicalInteractionLength(
							   const G4Track& aTrack,
							   G4double prevStepSize,
							   G4ForceCondition* cond)
{
  G4double me=electron_mass_c2/c_squared;
  G4double velLong=5324.2077*m/s;
  G4double mc, l0, ksound, kmax, velocity;
  if (aTrack.GetParticleDefinition()->GetParticleName() == "G4CMPDriftElectron")
  {
    G4ThreeVector T = G4ThreeVector(sqrt(.118/1.588), sqrt(.118/.081), sqrt(.118/.081));
    mc=.118*me;
    l0 = 257e-6*m;
    G4RotationMatrix trix;

    int valley = G4CMPValleyTrackMap::GetInstance()->GetValley(aTrack);
    
  switch(valley){
      case 1:
	  trix = G4RotationMatrix(G4Rep3x3( 1.0/sqrt(3.0),  1.0/sqrt(3.0),  1.0/sqrt(3.0), 
			  -1.0/sqrt(2.0),  1.0/sqrt(2.0),  0.0, 
			  -1.0/sqrt(6.0), -1.0/sqrt(6.0),  sqrt(2.0/3.0) ));
          break;
      case 2:
	  trix = G4RotationMatrix(G4Rep3x3(-1.0/sqrt(3.0),  1.0/sqrt(3.0),  1.0/sqrt(3.0), 
		 -1.0/sqrt(2.0), -1.0/sqrt(2.0),  0.0, 
		  1.0/sqrt(6.0), -1.0/sqrt(6.0),  sqrt(2.0/3.0) ));
          break;
      case 3:
	  trix = G4RotationMatrix(G4Rep3x3(-1.0/sqrt(3.0), -1.0/sqrt(3.0),  1.0/sqrt(3.0), 
		  1.0/sqrt(2.0), -1.0/sqrt(2.0),  0.0, 
		  1.0/sqrt(6.0),  1.0/sqrt(6.0),  sqrt(2.0/3.0) ));
          break;
      case 4:
	  trix = G4RotationMatrix(G4Rep3x3( 1.0/sqrt(3.0), -1.0/sqrt(3.0),  1.0/sqrt(3.0), 
		  1.0/sqrt(2.0),  1.0/sqrt(2.0),  0.0, 
		 -1.0/sqrt(6.0),  1.0/sqrt(6.0),  sqrt(2.0/3.0) ));
          break;
  }
    
    G4RotationMatrix mInv = 
		trix.inverse()*G4Rep3x3(1/1.588/me,   0.0    , 0.0,
							0.0     , 1/.081/me, 0.0, 
							0.0     ,   0.0    , 1/.081/me)
							*trix;
  
    valleyToNormal= G4AffineTransform(trix);
    normalToValley= G4AffineTransform(trix).Inverse();  
    G4FieldManager* fMan = 
	G4TransportationManager::GetTransportationManager()->GetFieldManager();
    const G4Field* field = fMan->GetDetectorField();

    G4ThreeVector posVec = aTrack.GetPosition();
    G4double position[4] = {posVec[0],posVec[1],posVec[2],0};
    G4double fieldVal[6];

    field->GetFieldValue(position,fieldVal);
    G4ThreeVector Efield = G4ThreeVector(fieldVal[3], fieldVal[4], fieldVal[5]);
    G4ThreeVector Efield_valley = normalToValley.TransformPoint(Efield);
    G4ThreeVector Efield_HV = G4ThreeVector( Efield_valley[0]*T[0], 
					     Efield_valley[1]*T[1], 
					     Efield_valley[2]*T[2]);
					     
    G4ThreeVector k = aTrack.GetMomentum()/hbarc;
    G4ThreeVector k_valley = normalToValley.TransformPoint(k);
    G4ThreeVector k_HV= G4ThreeVector( k_valley[0]*T[0],
					k_valley[1]*T[1], 
					k_valley[2]*T[2]);
    G4ThreeVector v_valley;
    v_valley[0] = hbar_Planck*k_valley[0]/1.588/me;
    v_valley[1] = hbar_Planck*k_valley[1]/.081/me;
    v_valley[2] = hbar_Planck*k_valley[2]/.081/me;
    G4ThreeVector v = valleyToNormal.TransformPoint(v_valley);
    
//     v_valley = k_HV*hbar_Planck/mc;
//     v_valley[0] *= T[0];
//     v_valley[1] *= T[1];
//     v_valley[2] *= T[2];
//     v = valleyToNormal.TransformPoint(v_valley);
    
    //v = k*hbar_Planck/mc;

    ksound = velLong*mc/hbar_Planck;
    //if (k_HV.mag()>ksound)
//	return DBL_MAX;
    kmax = 28.0*pow(m/volt, 
		    1.0/3.0)*ksound*pow(Efield_HV.mag()/10.0, 1.0/3.0);
    velocity = v.mag();
  }
  else
  {
    mc=.35*me;
    l0 = 108e-6*m;
    
    G4FieldManager* fMan = 
	G4TransportationManager::GetTransportationManager()->GetFieldManager();
    const G4Field* field = fMan->GetDetectorField();

    G4ThreeVector posVec = aTrack.GetPosition();
    G4double position[4] = {posVec[0],posVec[1],posVec[2],0};
    G4double fieldVal[6];

    field->GetFieldValue(position,fieldVal);
    G4ThreeVector Efield = G4ThreeVector(fieldVal[3], fieldVal[4], fieldVal[5]);

    ksound = velLong*mc/hbar_Planck;
    kmax = 14.72*pow(m/volt, 
		    1.0/3.0)*ksound*pow(Efield.mag()/10.0, 1.0/3.0);
    velocity = aTrack.GetStep()->GetPostStepPoint()->GetVelocity();
  }

    //set condition to "Forced"
    *cond = Forced;

    G4double dt =  1.0 / (
		       2 * velLong / (3*l0)
		       * (kmax / ksound) * (kmax / ksound)
		       * ((1- ksound/kmax))
		       * ((1- ksound/kmax))
		       * ((1- ksound/kmax))
		       );
    //G4cout << velocity*dt/m << G4endl;
    return velocity*dt;
}

G4VParticleChange* G4CMPTimeStepper::PostStepDoIt(
					  const G4Track& aTrack,
					  const G4Step& aStep
					  )
{
   //G4cout <<  "Time Step" <<  G4endl;
  aParticleChange.Initialize(aTrack);
  /*
  G4cout<<"\nG4CMPTimeStepper::PostStepDoIt: Delta time: " 
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

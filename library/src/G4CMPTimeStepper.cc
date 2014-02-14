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
  G4double velLong=5324.2077*m/s;
  G4double me = electron_mass_c2/c_squared;
  G4double mc_e=.118*me;
  G4double l0_e = 257e-6*m;
  G4ThreeVector T = G4ThreeVector(sqrt(.118/1.588), sqrt(.118/.081), sqrt(.118/.081));

//  valleyToNormal = G4AffineTransform(trix);
//  normalToValley = G4AffineTransform(trix).Inverse();

  G4FieldManager* fMan =
          G4TransportationManager::GetTransportationManager()->GetFieldManager();
  const G4Field* field = fMan->GetDetectorField();

  G4double position[4] = {0.0, 0.0, 0.0, 0.0};
  G4double fieldVal[6];

  field->GetFieldValue(position,fieldVal);
  G4ThreeVector Efield = G4ThreeVector(fieldVal[3], fieldVal[4], fieldVal[5]);
//  G4ThreeVector Efield_valley = normalToValley.TransformPoint(Efield);
//  G4ThreeVector Efield_HV = G4ThreeVector( Efield_valley[0]*T[0],
//                                           Efield_valley[1]*T[1],
//                                           Efield_valley[2]*T[2]);

  G4double ksound_e = velLong*mc_e/hbar_Planck;

  G4double kmaxElec = 28.0 * pow(m/volt, 1.0/3.0)
                      * ksound_e * pow(Efield.mag(), 1.0/3.0);

  dt_e =  1.0 / ( 2 * velLong / (3*l0_e)
           * (kmaxElec / ksound_e) * (kmaxElec / ksound_e)
           * ((1- ksound_e/kmaxElec))
           * ((1- ksound_e/kmaxElec))
           * ((1- ksound_e/kmaxElec)) );


  G4double mc_h=.35*me;
  G4double l0_h = 108e-6*m;

  G4double ksound_h = velLong*mc_h/hbar_Planck;

  G4double kmaxHole = 14.72 * pow(m/volt, 1.0/3.0)
                      * ksound_h * pow(Efield.mag()/10.0, 1.0/3.0);

  dt_h =  1.0 / ( 2 * velLong / (3*l0_h)
           * (kmaxHole / ksound_h) * (kmaxHole / ksound_h)
           * ((1- ksound_h/kmaxHole))
           * ((1- ksound_h/kmaxHole))
           * ((1- ksound_h/kmaxHole)) );

  if (verboseLevel>0) {
    G4cout << GetProcessName() << " is created " << G4endl;
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4CMPTimeStepper::~G4CMPTimeStepper()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CMPTimeStepper::G4CMPTimeStepper(G4CMPTimeStepper& right)
  : G4VProcess(right), dt_e(right.dt_e), dt_h(right.dt_h)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4CMPTimeStepper::
PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
				     G4double /*prevStepSize*/,
				     G4ForceCondition* /*cond*/) {
  // Only drifting electrons have special treatment
  if (aTrack.GetParticleDefinition() != G4CMPDriftElectron::Definition()) {
    G4double v = aTrack.GetStep()->GetPostStepPoint()->GetVelocity();
    return v*dt_h;
  }

  G4double me = electron_mass_c2/c_squared;
  G4RotationMatrix trix;

  int valley = G4CMPValleyTrackMap::GetInstance()->GetValley(aTrack);
  switch (valley) { 
  case 1:
    trix = G4RotationMatrix(
		G4Rep3x3( 1.0/sqrt(3.0),  1.0/sqrt(3.0),  1.0/sqrt(3.0),
			 -1.0/sqrt(2.0),  1.0/sqrt(2.0),  0.0,
			 -1.0/sqrt(6.0), -1.0/sqrt(6.0),  sqrt(2.0/3.0) ));
    break;
  case 2:
    trix = G4RotationMatrix(
		G4Rep3x3(-1.0/sqrt(3.0),  1.0/sqrt(3.0),  1.0/sqrt(3.0),
			 -1.0/sqrt(2.0), -1.0/sqrt(2.0),  0.0,
			  1.0/sqrt(6.0), -1.0/sqrt(6.0),  sqrt(2.0/3.0) ));
    break;
  case 3:
    trix = G4RotationMatrix(
		G4Rep3x3(-1.0/sqrt(3.0), -1.0/sqrt(3.0),  1.0/sqrt(3.0),
			  1.0/sqrt(2.0), -1.0/sqrt(2.0),  0.0,
			  1.0/sqrt(6.0),  1.0/sqrt(6.0),  sqrt(2.0/3.0) ));
    break;
  case 4:
    trix = G4RotationMatrix(
		G4Rep3x3( 1.0/sqrt(3.0), -1.0/sqrt(3.0),  1.0/sqrt(3.0),
			  1.0/sqrt(2.0),  1.0/sqrt(2.0),  0.0,
			 -1.0/sqrt(6.0),  1.0/sqrt(6.0),  sqrt(2.0/3.0) ));
    break;
  }

  G4RotationMatrix mInv = trix.inverse()
    * G4Rep3x3(1/1.588/me,   0.0    , 0.0,
	         0.0     , 1/.081/me, 0.0,
	         0.0     ,   0.0    , 1/.081/me)
    * trix;

  G4ThreeVector k = aTrack.GetMomentum()/hbarc;
  G4ThreeVector k_valley = normalToValley.TransformPoint(k);
  
  G4ThreeVector v_valley = hbar_Planck*(mInv*k_valley);
  G4ThreeVector v = valleyToNormal.TransformPoint(v_valley);
  
  return v.mag()*dt_e;
}


G4VParticleChange* G4CMPTimeStepper::PostStepDoIt(const G4Track& aTrack,
						  const G4Step& /*aStep*/) {
  aParticleChange.Initialize(aTrack);
  return &aParticleChange;
}

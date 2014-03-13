#include "G4CMPInterValleyScattering.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPValleyTrackMap.hh"
#include "G4Field.hh"
#include "G4FieldManager.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalVolume.hh"
#include "G4PhysicalConstants.hh"
#include "G4RandomDirection.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4VParticleChange.hh"
#include "G4VPhysicalVolume.hh"
#include "Randomize.hh"

#include "math.h"

G4CMPInterValleyScattering::G4CMPInterValleyScattering()
  : G4VDiscreteProcess("InterValleyScattering"), G4CMPProcessUtils() {
  //E_0 from Edelweiss experiment LTD 14. This value is tuned for
  //EDELWEISS/CDMS crystals with small (~V/m) electric fields. It 
  //may still be fine for larger electric fields (~10-100V/m) 
  //however there is no experimental data as of now. 
  E_0_ED_203 = 217.0 ; // V/m
    
  if (verboseLevel>1) {
    G4cout<<GetProcessName()<<" is created "<<G4endl;
  }
}

G4CMPInterValleyScattering::~G4CMPInterValleyScattering()
{ ; }

G4CMPInterValleyScattering::G4CMPInterValleyScattering(G4CMPInterValleyScattering& right)
: G4VDiscreteProcess(right)
{ ; }

G4double 
G4CMPInterValleyScattering::GetMeanFreePath(const G4Track& aTrack,
					    G4double,
					    G4ForceCondition* condition) {
  G4StepPoint* stepPoint  = aTrack.GetStep()->GetPostStepPoint();
  G4double velocity = stepPoint->GetVelocity();
  *condition = NotForced;
  
  const G4RotationMatrix& trix = GetValley(aTrack);
  normalToValley = G4AffineTransform(trix);
  valleyToNormal = G4AffineTransform(trix).Inverse();
  
  //G4VPhysicalVolume* pVol = aTrack.GetVolume();
  //G4LogicalVolume* lVol = pVol->GetLogicalVolume();
  
  G4FieldManager* fMan = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  
  //If there is no field, there is no IV scattering... but then there
  //is no e-h transport either...
  if (!fMan->DoesFieldExist()) return DBL_MAX;
  
  //Getting the field.
  const G4Field* field = fMan->GetDetectorField();
  
  G4ThreeVector pos = aTrack.GetPosition();
  G4double  posVec[4] = { pos.x(), pos.y(), pos.z(), 0. };
  
  G4double fieldValue[6];
  
  field->GetFieldValue(posVec,fieldValue);
  
  G4ThreeVector fieldVector = G4ThreeVector(fieldValue[3]/(volt/m),
					    fieldValue[4]/(volt/m), 
					    fieldValue[5]/(volt/m));
  
  G4ThreeVector fieldVectorLLT = 
    normalToValley.TransformAxis(fieldVector).unit();
  G4ThreeVector fieldDirHV = G4ThreeVector(fieldVectorLLT.getX()*(1/1.2172),
					   fieldVectorLLT.getY()*(1/1.2172),
					   fieldVectorLLT.getZ()*(1/0.27559)
					   );
  
  G4double fieldMagHV = fieldDirHV.mag();
  
  //setting the E-field parameter as outline in EDELSWEISS LTD-14.
  E_0_ED_201 = fieldMagHV; 
  
  G4double mfp = velocity * 6.72e-6 * 
    pow((E_0_ED_203 * E_0_ED_203 + abs(E_0_ED_201) * abs(E_0_ED_201)), 3.24/2.0 );
  
  return mfp;
}

G4VParticleChange* 
G4CMPInterValleyScattering::PostStepDoIt(const G4Track& aTrack, 
					 const G4Step& /*aStep*/) {
  // picking a new valley at random if IV-scattering process was triggered
  int valley = ChooseValley();
    
  // Assigning a new valley...
  G4CMPValleyTrackMap::GetInstance()->SetValley(aTrack, valley);
  
  // No change to actual momentum, just following a different valley
  aParticleChange.Initialize(aTrack);  
  
  ResetNumberOfInteractionLengthLeft();    
  return &aParticleChange;
}

G4bool G4CMPInterValleyScattering::IsApplicable(const G4ParticleDefinition& aPD)
{
  return (&aPD==G4CMPDriftElectron::Definition());
}

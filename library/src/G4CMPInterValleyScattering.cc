// $Id$
//
// 20140324  Drop hard-coded IV scattering parameters; get from lattice
// 20140324  Restore Z-axis mass tensor
// 20140331  Add required process subtype code
// 20140418  Drop local valley transforms, use lattice functions instead
// 20140429  Recompute kinematics relative to new valley
// 20140908  Allow IV scatter to change momentum by conserving energy
// 20150109  Revert IV scattering to preserve momentum
// 20150112  Follow renaming of "SetNewKinematics" to FillParticleChange
// 20150122  Use verboseLevel instead of compiler flag for debugging

#include "G4CMPInterValleyScattering.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPValleyTrackMap.hh"
#include "G4CMPFieldManager.hh"
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
#include "G4VParticleChange.hh"
#include "G4VPhysicalVolume.hh"
#include "Randomize.hh"
#include "math.h"

G4CMPInterValleyScattering::G4CMPInterValleyScattering()
  : G4CMPVDriftProcess("G4CMPInterValleyScattering", fInterValleyScattering) {;}

G4CMPInterValleyScattering::~G4CMPInterValleyScattering() {;}


G4double 
G4CMPInterValleyScattering::GetMeanFreePath(const G4Track& aTrack,
					    G4double,
					    G4ForceCondition* condition) {
  *condition = NotForced;

  // Get electric field associated with current volume, if any
  G4FieldManager* fMan =
    aTrack.GetVolume()->GetLogicalVolume()->GetFieldManager();
  
  //If there is no field, there is no IV scattering... but then there
  //is no e-h transport either...
  if (!fMan || !fMan->DoesFieldExist()) return DBL_MAX;

  G4ThreeVector p_local = GetLocalMomentum(aTrack);

  G4int valley = GetValleyIndex(aTrack);
  G4double velocity = theLattice->MapPtoV_el(valley, p_local).mag();
  
  G4double posVec[4] = { 4*0. };
  GetLocalPosition(aTrack, posVec);

  const G4Field* field = fMan->GetDetectorField();
  G4double fieldValue[6];
  field->GetFieldValue(posVec,fieldValue);

  G4ThreeVector fieldVector(fieldValue[3], fieldValue[4], fieldValue[5]);

  if (verboseLevel > 1) {
    G4cout << "IV local position (" << posVec[0] << "," << posVec[1] << ","
	   << posVec[2] << ")\n field " << fieldVector/volt*cm << " V/cm"
	   << "\n magnitude " << fieldVector.mag()/volt*cm << " V/cm toward "
	   << fieldVector.cosTheta() << " z" << G4endl;
  }

  // Find E-field in HV space: rotate into valley, then apply HV tansform.
  // NOTE:  Separate steps to avoid matrix-matrix multiplications
  fieldVector *= GetValley(aTrack);
  fieldVector *= theLattice->GetSqrtInvTensor();
  fieldVector /= volt/m;			// Strip units for MFP below

  if (verboseLevel > 1) {
    G4cout << " in HV space " << fieldVector*0.01 << " ("
	   << fieldVector.mag()*0.01 << ") V/cm" << G4endl;
  }

  // Compute mean free path per Edelweiss LTD-14 paper
  G4double E_0 = theLattice->GetIVField();
  G4double mfp = velocity / ( theLattice->GetIVRate() *
    pow((E_0*E_0 + fieldVector.mag2()), theLattice->GetIVExponent()/2.0) );

  if (verboseLevel > 1) G4cout << "IV MFP = " << mfp/m << G4endl;
  return mfp;
}

G4VParticleChange* 
G4CMPInterValleyScattering::PostStepDoIt(const G4Track& aTrack, 
					 const G4Step& aStep) {
  aParticleChange.Initialize(aTrack); 
  G4StepPoint* postStepPoint = aStep.GetPostStepPoint();
  
  // Do nothing when step limit reached or leaving volume
  if (verboseLevel > 0) {
    G4cout << GetProcessName() << "::PostStepDoIt: Step limited by process "
	   << postStepPoint->GetProcessDefinedStep()->GetProcessName()
	   << G4endl;
  }

  if (postStepPoint->GetStepStatus()==fGeomBoundary) {
    return &aParticleChange;
  }

  // Get track's energy in current valley
  G4ThreeVector p = GetLocalMomentum(aTrack);
  G4int valley = GetValleyIndex(aTrack);
  p = theLattice->MapPtoK_valley(valley, p); // p is actually k now

  // picking a new valley at random if IV-scattering process was triggered
  valley = ChooseValley();
  trackVmap->SetValley(aTrack, valley);

  p = theLattice->MapK_valleyToP(valley, p); // p is p again
  RotateToGlobalDirection(p);

  // Adjust track kinematics for new valley
  FillParticleChange(valley, p);

  ResetNumberOfInteractionLengthLeft();    
  return &aParticleChange;
}

G4bool G4CMPInterValleyScattering::IsApplicable(const G4ParticleDefinition& aPD)
{
  return (&aPD==G4CMPDriftElectron::Definition());
}

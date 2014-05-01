// $Id$
//
// 20140324  Drop hard-coded IV scattering parameters; get from lattice
// 20140324  Restore Z-axis mass tensor
// 20140331  Add required process subtype code
// 20140418  Drop local valley transforms, use lattice functions instead
// 20140429  Recompute kinematics relative to new valley

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
#include "G4VParticleChange.hh"
#include "G4VPhysicalVolume.hh"
#include "Randomize.hh"
#include "math.h"

G4CMPInterValleyScattering::G4CMPInterValleyScattering()
  : G4CMPVDriftProcess("InterValleyScattering", fInterValleyScattering) {;}

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
  if (!fMan->DoesFieldExist()) return DBL_MAX;
  
  G4StepPoint* stepPoint  = aTrack.GetStep()->GetPostStepPoint();
  G4double velocity = stepPoint->GetVelocity();
  
  G4double posVec[4] = { 4*0. };
  GetLocalPosition(aTrack, posVec);

  const G4Field* field = fMan->GetDetectorField();
  G4double fieldValue[6];
  field->GetFieldValue(posVec,fieldValue);

  G4ThreeVector fieldVector(fieldValue[3], fieldValue[4], fieldValue[5]);

  // Find E-field in HV space by rotating into valley and then applying HV tansform.
  // Also have to strip Efield units for use in MFP calculation.
  fieldVector = theLattice->GetSqrtInvTensor() * GetValley(aTrack) * fieldVector/volt*m;

  // Compute mean free path per Edelweiss LTD-14 paper
  G4double E_0 = theLattice->GetIVField();
  G4double mfp = velocity * theLattice->GetIVRate()/s *
    pow((E_0*E_0 + fieldVector.mag2()), theLattice->GetIVExponent()/2.0);

#ifdef G4CMP_DEBUG
  G4cout << "IV MFP = " << mfp/m << G4endl;
#endif
  return mfp;
}

G4VParticleChange* 
G4CMPInterValleyScattering::PostStepDoIt(const G4Track& aTrack, 
					 const G4Step& /*aStep*/) {
  // picking a new valley at random if IV-scattering process was triggered
  int valley = ChooseValley();
    
  // Assigning a new valley...
  trackVmap->SetValley(aTrack, valley);
  
  // Adjust kinetic energy and "effective mass" for new valley frame
  // NOTE:  These _should_ be the same if scattering is conservative
  aParticleChange.Initialize(aTrack);  
  SetNewKinematics(valley, aTrack.GetMomentum());

  ResetNumberOfInteractionLengthLeft();    
  return &aParticleChange;
}

G4bool G4CMPInterValleyScattering::IsApplicable(const G4ParticleDefinition& aPD)
{
  return (&aPD==G4CMPDriftElectron::Definition());
}

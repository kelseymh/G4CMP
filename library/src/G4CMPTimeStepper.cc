// $Id$
//
// 20140313  Introduce multiple inheritance from G4CMPProcessUtils, also
//	     use proper TransformAxis() on vectors, *not* TransformPoint()
//	     Add wrapper function to compute individual time steps
// 20140314  Fetch propagation parameters from lattice, instead of hardwired
// 20140324  Migrate to use of volume-local field, do coordinate transforms
// 20140325  Move most of time-step calculation to G4CMPProcessUtils

#include "G4CMPTimeStepper.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4CMPValleyTrackMap.hh"
#include "G4Field.hh"
#include "G4FieldManager.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalVolume.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4UserLimits.hh"
#include "G4VParticleChange.hh"
#include "G4VPhysicalVolume.hh"
#include <fstream>
#include <iostream>
#include <math.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CMPTimeStepper::G4CMPTimeStepper()
  : G4CMPVDriftProcess("G4CMPTimeStepper"), dt_e(0.), dt_h(0.) {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CMPTimeStepper::~G4CMPTimeStepper() {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4CMPTimeStepper::IsApplicable(const G4ParticleDefinition &pd) {
  return (&pd == G4CMPDriftElectron::Definition() ||
	  &pd == G4CMPDriftHole::Definition());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4CMPTimeStepper::
PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
				     G4double /*prevStepSize*/,
				     G4ForceCondition* /*cond*/) {
  ComputeTimeSteps(aTrack);

  // Only drifting electrons have special treatment
  if (aTrack.GetParticleDefinition() != G4CMPDriftElectron::Definition()) {
    G4double v = aTrack.GetStep()->GetPostStepPoint()->GetVelocity();
#ifdef G4CMP_DEBUG
    G4cout << "TS hole = " << (v*dt_h)/m << G4endl;
#endif
    return v*dt_h;
  }

  const G4RotationMatrix& trix = GetValley(aTrack);
  valleyToNormal = G4AffineTransform(trix);
  normalToValley = G4AffineTransform(trix).Inverse();

  G4ThreeVector k = aTrack.GetMomentum()/hbarc;
  G4ThreeVector k_valley = normalToValley.TransformAxis(k);
  
  G4ThreeVector v_valley = hbar_Planck * (theLattice->GetMInvTensor()*k_valley);
  G4ThreeVector v = valleyToNormal.TransformAxis(v_valley);
  
#ifdef G4CMP_DEBUG
  G4cout << "TS electron = " << (v.mag()*dt_e)/m << G4endl;
#endif  
  return v.mag()*dt_e;
}


G4VParticleChange* G4CMPTimeStepper::PostStepDoIt(const G4Track& aTrack,
						  const G4Step& /*aStep*/) {
  aParticleChange.Initialize(aTrack);
  return &aParticleChange;
}


// Compute dt_e, dt_h and valley rotations at current location

void G4CMPTimeStepper::ComputeTimeSteps(const G4Track& aTrack) {
  G4FieldManager* fMan =
    aTrack.GetVolume()->GetLogicalVolume()->GetFieldManager();

  if (!fMan || !fMan->DoesFieldExist()) {	// No field, no special action
    dt_e = 3.*l0_e/velLong;
    dt_h = 3.*l0_h/velLong;
    return;
  }

  G4double position[4] = {4*0.};
  GetLocalPosition(aTrack, position);

  const G4Field* field = fMan->GetDetectorField();
  G4double fieldVal[6];
  field->GetFieldValue(position,fieldVal);

  G4ThreeVector Efield(fieldVal[3], fieldVal[4], fieldVal[5]);

  dt_e = TimeStepInField(Efield.mag()/10., 28.0, l0_e);
  dt_h = TimeStepInField(Efield.mag()/10., 14.72, l0_h);
}

// Compute time step for electrons or holes, pre-simplified expression

G4double G4CMPTimeStepper::TimeStepInField(G4double Emag, G4double coeff, G4double l0) const {
  // Maximum "Mach number" (carrier speed vs. sound) at specified field
  G4double maxMach = coeff * pow(Emag*(m/volt), 1./3.);

  return (0.5 * ChargeCarrierTimeStep(maxMach, l0));
}

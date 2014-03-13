// $Id$
//
// 20140313  Introduce multiple inheritance from G4CMPProcessUtils, also
//	     use proper TransformAxis() on vectors, *not* TransformPoint()

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
#include "G4TransportationManager.hh"
#include "G4UserLimits.hh"
#include "G4VParticleChange.hh"
#include "G4VPhysicalVolume.hh"
#include <fstream>
#include <iostream>
#include <math.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CMPTimeStepper::G4CMPTimeStepper()
  : G4VDiscreteProcess("G4CMPTimeStepper", fGeneral), G4CMPProcessUtils(),
    dt_e(0.), dt_h(0.),
    velLong(5324.2077*m/s), me(electron_mass_c2/c_squared),
    mc_e(.118*me), l0_e(257e-6*m), ksound_e(velLong*mc_e/hbar_Planck),
    mc_h(.350*me), l0_h(108e-6*m), ksound_h(velLong*mc_h/hbar_Planck) {
  if (verboseLevel>0) {
    G4cout << GetProcessName() << " is created " << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CMPTimeStepper::~G4CMPTimeStepper()
{;}

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
    return v*dt_h;
  }

  const G4RotationMatrix& trix = GetValley(aTrack);
  valleyToNormal = G4AffineTransform(trix);
  normalToValley = G4AffineTransform(trix).Inverse();

  G4RotationMatrix mInv = trix.inverse()*theLattice->GetElectronMass()*trix;

  G4ThreeVector k = aTrack.GetMomentum()/hbarc;
  G4ThreeVector k_valley = normalToValley.TransformAxis(k);
  
  G4ThreeVector v_valley = hbar_Planck*(mInv*k_valley);
  G4ThreeVector v = valleyToNormal.TransformAxis(v_valley);
  
  return v.mag()*dt_e;
}


G4VParticleChange* G4CMPTimeStepper::PostStepDoIt(const G4Track& aTrack,
						  const G4Step& /*aStep*/) {
  aParticleChange.Initialize(aTrack);
  return &aParticleChange;
}


// Compute dt_e, dt_h and valley rotations at current location

void G4CMPTimeStepper::ComputeTimeSteps(const G4Track& aTrack) {
  G4ThreeVector pos = aTrack.GetPosition();
  G4double position[4] = {pos.x(), pos.y(), pos.z(), 0.0};

  G4FieldManager* fMan =
    G4TransportationManager::GetTransportationManager()->GetFieldManager();
  if (!fMan->DoesFieldExist()) {
    dt_e = 3.*l0_e/velLong;
    dt_h = 3.*l0_h/velLong;
    return;
  }

  const G4Field* field = fMan->GetDetectorField();

  G4double fieldVal[6];

  field->GetFieldValue(position,fieldVal);
  G4ThreeVector Efield = G4ThreeVector(fieldVal[3], fieldVal[4], fieldVal[5]);

  G4double kmaxElec = 28.0 * ksound_e * pow(Efield.mag()*(m/volt), 1.0/3.0);

  dt_e =  1.0 / ( 2 * velLong / (3*l0_e)
           * (kmaxElec / ksound_e) * (kmaxElec / ksound_e)
           * ((1- ksound_e/kmaxElec))
           * ((1- ksound_e/kmaxElec))
           * ((1- ksound_e/kmaxElec)) );

  G4double kmaxHole = 14.72 * ksound_h * pow(Efield.mag()*(m/volt)/10.0, 1.0/3.0);

  dt_h =  1.0 / ( 2 * velLong / (3*l0_h)
           * (kmaxHole / ksound_h) * (kmaxHole / ksound_h)
           * ((1- ksound_h/kmaxHole))
           * ((1- ksound_h/kmaxHole))
           * ((1- ksound_h/kmaxHole)) );
}

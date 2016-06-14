/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
//
// 20140313  Introduce multiple inheritance from G4CMPProcessUtils, also
//	     use proper TransformAxis() on vectors, *not* TransformPoint()
//	     Add wrapper function to compute individual time steps
// 20140314  Fetch propagation parameters from lattice, instead of hardwired
// 20140324  Migrate to use of volume-local field, do coordinate transforms
// 20140325  Move most of time-step calculation to G4CMPProcessUtils
// 20140331  Add required process subtype code
// 20140418  Remove explicit valley transforms, use lattice function
// 20140429  Adjust "effective mass" and energy based on end-of-step momentum
// 20150112  Follow renaming of "SetNewKinematics" to FillParticleChange, drop
//	     redundant IsApplicable()
// 20150122  Use verboseLevel instead of compiler flag for debugging

#include "G4CMPTimeStepper.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4CMPTrackInformation.hh"
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
#include <math.h>

G4CMPTimeStepper::G4CMPTimeStepper()
  : G4CMPVDriftProcess("G4CMPTimeStepper", fTimeStepper) {;}


G4CMPTimeStepper::~G4CMPTimeStepper() {;}


G4double G4CMPTimeStepper::
PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
				     G4double /*prevStepSize*/,
				     G4ForceCondition* /*cond*/) {
  G4CMPTrackInformation* trackInfo = static_cast<G4CMPTrackInformation*>(
    aTrack.GetAuxiliaryTrackInformation(fPhysicsModelID));
  G4ThreeVector x0 = GetLocalPosition(aTrack);
  G4ThreeVector v0 = GetLocalVelocityVector(aTrack);
  G4ThreeVector k0 = GetLocalWaveVector(aTrack);
  G4ThreeVector Efield0 = GetLocalEffectiveEField(aTrack, x0);
  G4double mass = trackInfo->GetEffectiveMass()*c_squared;
  G4double q = aTrack.GetDynamicParticle()->GetCharge()*eplus;
  G4double kSound = theLattice->GetSoundSpeed();
  G4double l0 = trackInfo->GetScatterLength();

  // Leman-Verlet method:
  G4double dt0 = ComputeTimeSteps(aTrack);
  G4ThreeVector x1 = x0 + v0*dt0 + 0.5*dt0*dt0*q*Efield0/mass;
  G4ThreeVector Efield1 = GetLocalEffectiveEField(aTrack, x1);
  G4ThreeVector k1 = k0 + 0.5*q/hbarc*dt0*(Efield0 + Efield1);
  if (k0.mag() <= kSound || k1.mag() <= kSound) {
    return v0.mag()*dt0;
  }

  G4double tau0 = ChargeCarrierTimeStep(k0.mag()/kSound, l0);
  G4double tau1 = ChargeCarrierTimeStep(k1.mag()/kSound, l0);
  G4double a0 = 1./tau0;
  G4double a1 = 1./dt0*(1./tau1 - 1./tau0);
  G4double arg = a0*a0 - 2.*a1*std::log(G4UniformRand());
  if (arg < 0. || a1 < DBL_MIN)
    return v0.mag()*dt0;
  G4double dt1 = (std::sqrt(arg) - a0)/a1;
  if (dt1 < DBL_MIN) dt1 = DBL_MIN;

  if (dt1 < dt0) {
    x1 = x0 + v0*dt1 + 0.5*dt1*dt1*q/mass*Efield0;
    Efield1 = GetLocalEffectiveEField(aTrack, x1);
    k1 = k0 + 0.5*q/hbarc*dt1*(Efield0 + Efield1);
    if (k1.mag() <= kSound) {
      return v0.mag()*dt0;
    }
    tau1 = ChargeCarrierTimeStep(k1.mag()/kSound, l0);
    a1 = 1./dt1*(1./tau1 - 1./tau0);
    arg = a0*a0 - 2.*a1*std::log(G4UniformRand());
    if (arg < 0. || a1 < DBL_MIN)
      return v0.mag()*dt0;
    dt1 = (std::sqrt(arg) - a0)/a1;
    if (dt1 < DBL_MIN) dt1 = DBL_MIN;
  }
   G4double gpil = 0.5*((v0+ConvertWaveVectorToVelocityVector(aTrack, k1))*dt1).mag();

  if (verboseLevel > 1) {
      G4cout << "TimeStepper " << aTrack.GetDefinition()->GetParticleName() <<
                " GPIL = " << gpil/m << G4endl;
  }

  return gpil;
}


G4VParticleChange* G4CMPTimeStepper::PostStepDoIt(const G4Track& aTrack,
						  const G4Step& aStep) {
  aParticleChange.Initialize(aTrack);

  // Adjust mass and kinetic energy using end-of-step momentum
  G4ThreeVector pfinal = GetGlobalMomentum(aTrack);
  FillParticleChange(GetValleyIndex(aTrack), pfinal);

  return &aParticleChange;
//  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

// Compute dt_e, dt_h and valley rotations at current location

G4double G4CMPTimeStepper::ComputeTimeSteps(const G4Track& aTrack) {
  G4CMPTrackInformation* trackInfo = static_cast<G4CMPTrackInformation*>(
    aTrack.GetAuxiliaryTrackInformation(fPhysicsModelID));
  G4double l0 = trackInfo->GetScatterLength();

  G4ThreeVector Efield = GetLocalEffectiveEField(aTrack);
  if (Efield.mag() <= DBL_MIN) // No field, no special action
    return 3.*l0/velLong;

  G4double timeStepParam = 0.;
  if (aTrack.GetParticleDefinition() == G4CMPDriftElectron::Definition()) {
    timeStepParam = 28.0;
  } else if (aTrack.GetParticleDefinition() == G4CMPDriftHole::Definition()) {
    timeStepParam = 14.72;
  } else {
    // TODO: Write an exception
  }

  return TimeStepInField(Efield.mag()/10., timeStepParam, l0);
}

// Compute time step for electrons or holes, pre-simplified expression

G4double G4CMPTimeStepper::TimeStepInField(G4double Emag, G4double coeff, G4double l0) const {
  // Maximum "Mach number" (carrier speed vs. sound) at specified field
  G4double maxMach = coeff * pow(Emag*(m/volt), 1./3.);

  return (0.5 * ChargeCarrierTimeStep(maxMach, l0));
}

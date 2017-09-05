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
// 20160624  Use GetTrackInfo() accessor
// 20161114  Use new G4CMPDriftTrackInfo
// 20170602  Use G4CMPUtils for track identity functions
// 20170806  Swap GPIL and MFP functions to work with G4CMPVProcess base

#include "G4CMPTimeStepper.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4CMPDriftTrackInfo.hh"
#include "G4CMPTrackUtils.hh"
#include "G4CMPUtils.hh"
#include "G4CMPVProcess.hh"
#include "G4CMPVScatteringRate.hh"
#include "G4Field.hh"
#include "G4FieldManager.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalVolume.hh"
#include "G4PhysicalConstants.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4UserLimits.hh"
#include "G4VParticleChange.hh"
#include "G4VPhysicalVolume.hh"
#include <math.h>

G4CMPTimeStepper::G4CMPTimeStepper()
  : G4CMPVDriftProcess("G4CMPTimeStepper", fTimeStepper) {;}

G4CMPTimeStepper::~G4CMPTimeStepper() {;}


// Get scattering rates from current track's processes

void G4CMPTimeStepper::LoadDataForTrack(const G4Track* aTrack) {
  G4CMPProcessUtils::LoadDataForTrack(aTrack);	// Common configuration

  lukeRate = ivRate = nullptr;			// Discard previous versions

  // Pointers can't be null since track has at least this process!
  const G4ProcessVector* pvec =
    aTrack->GetDefinition()->GetProcessManager()->GetPostStepProcessVector();

  if (verboseLevel>2)
    G4cout << "TimeStepper scanning " << pvec->size() << " processes"
	   << " for " << aTrack->GetDefinition()->GetParticleName() << G4endl;

  for (G4int i=0; i<pvec->size(); i++) {
    const G4CMPVProcess* cmpProc = dynamic_cast<G4CMPVProcess*>((*pvec)[i]);
    if (!cmpProc) continue;

    const G4String& pname = cmpProc->GetProcessName();
    if (verboseLevel>2) G4cout << pname << G4endl;

    if (pname == "G4CMPLukeScattering")      lukeRate = cmpProc->GetRateModel();
    if (pname == "G4CMPInterValleyScattering") ivRate = cmpProc->GetRateModel();
  }

  if (verboseLevel>1) {
    G4cout << "TimeStepper Found" << (lukeRate?" lukeRate":"")
	   << (ivRate?" ivRate":"") << G4endl;
  }
}


// Compute fixed "minimum distance" to avoid accelerating past Luke or IV

G4double G4CMPTimeStepper::GetMeanFreePath(const G4Track& aTrack, G4double,
					   G4ForceCondition* cond) {
  *cond = NotForced;

  G4double dt = ComputeTimeSteps(aTrack);
  G4double v = GetVelocity(aTrack);

  if (verboseLevel > 1) {
    G4cout << "TS " << (IsElectron()?"elec":"hole") << " = " << (v*dt)/m
	   << " m" << G4endl;
  }

  return v*dt;
}


// At end of step, recompute kinematics; important for electrons

G4VParticleChange* G4CMPTimeStepper::PostStepDoIt(const G4Track& aTrack,
						  const G4Step& aStep) {
  aParticleChange.Initialize(aTrack);

  // Adjust mass and kinetic energy using end-of-step momentum
  G4ThreeVector pfinal = GetGlobalMomentum(aTrack);
  FillParticleChange(GetValleyIndex(aTrack), pfinal);

  ClearNumberOfInteractionLengthLeft();		// All processes must do this!
  return &aParticleChange;
}


// Compute dt_e, dt_h and valley rotations at current location

G4double G4CMPTimeStepper::ComputeTimeSteps(const G4Track& aTrack) {
  G4double timeStepParam = 0.; G4double l0 = 0.;
  if (G4CMP::IsElectron(aTrack)) {
    timeStepParam = 28.0;
    l0 = theLattice->GetElectronScatter();
  } else if (G4CMP::IsHole(&aTrack)) {
    timeStepParam = 14.72;
    l0 = theLattice->GetHoleScatter();
  } else {
    // TODO: Write an exception
  }

  G4FieldManager* fMan =
    aTrack.GetVolume()->GetLogicalVolume()->GetFieldManager();
  if (!fMan || !fMan->DoesFieldExist()) {	// No field, no special action
    return ChargeCarrierTimeStep(0., l0);
  }

  G4double position[4] = { 4*0. };
  GetLocalPosition(aTrack, position);

  const G4Field* field = fMan->GetDetectorField();
  G4double fieldVal[6];
  field->GetFieldValue(position,fieldVal);

  G4ThreeVector Efield(fieldVal[3], fieldVal[4], fieldVal[5]);
  return TimeStepInField(Efield.mag()/10., timeStepParam, l0);
}


// Compute time step for electrons or holes, pre-simplified expression

G4double G4CMPTimeStepper::TimeStepInField(G4double Emag, G4double coeff, G4double l0) const {
  // Maximum "Mach number" (carrier speed vs. sound) at specified field
  G4double maxMach = coeff * pow(Emag*(m/volt), 1./3.);

  return (0.5 * ChargeCarrierTimeStep(maxMach, l0));
}

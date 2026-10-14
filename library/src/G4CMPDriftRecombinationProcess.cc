/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPDriftRecombinationProcess.cc
/// \brief Charge (electron-hole) recombination
//
// 20170620  M. Kelsey -- Follow interface changes in G4CMPSecondaryUtils
// 20170802  M. Kelsey -- Replace phonon production with G4CMPEnergyPartition
// 20180827  M. Kelsey -- Prevent partitioner from recomputing sampling factors
// 20210328  Modify above; compute direct-phonon sampling factor here
// 20250929  M. Kelsey -- Include residual kinetic energy in phonon release
// 20260514  G4CMP-517 -- Register NTL (Luke) rate model to use in computing
//	       threshold and energy-gain estimation for killing in flight.
// 20260903  G4CMP-517 -- Use new "recombinationScale" to modulate setting
//	       limits for "free flight" -- DBL_MAX restores old behaviour.

#include "G4CMPDriftRecombinationProcess.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4CMPEnergyPartition.hh"
#include "G4CMPFieldUtils.hh"
#include "G4CMPLukeEmissionRate.hh"
#include "G4CMPSecondaryUtils.hh"
#include "G4CMPSolidUtils.hh"
#include "G4CMPUtils.hh"
#include "G4LatticePhysical.hh"
#include "G4RandomDirection.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"
#include <vector>


// Constructor and destructor

G4CMPDriftRecombinationProcess::
G4CMPDriftRecombinationProcess(const G4String &name, G4CMPProcessSubType type)
  : G4CMPVDriftProcess(name, type), partitioner(new G4CMPEnergyPartition) {
  partitioner->UseDownsampling(false);		// Apply preset scaling factors
  UseRateModel(new G4CMPLukeEmissionRate);	// For threshold estimation
}

G4CMPDriftRecombinationProcess::~G4CMPDriftRecombinationProcess() {
  delete partitioner;
}


// Process actions

G4double 
G4CMPDriftRecombinationProcess::GetMeanFreePath(const G4Track& aTrack, G4double,
						G4ForceCondition* cond) {
  UpdateMeanFreePathForLatticeChangeover(aTrack);

  *cond = Forced;
  return DBL_MAX;
}

G4VParticleChange* 
G4CMPDriftRecombinationProcess::PostStepDoIt(const G4Track& aTrack,
					     const G4Step& aStep) {
  InitializeParticleChange(aTrack);

  if (!ReadyToRecombine(aTrack)) {	// Evaluate _after_ Transportation
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }

  if (verboseLevel) {
    G4cout << GetProcessName() << "::PostStepDoIt: " << G4endl
           << aTrack.GetDefinition()->GetParticleName() << " "
	   << aTrack.GetKineticEnergy()/eV << " eV "
	   << "reabsorbed by lattice @ " << aTrack.GetPosition() << G4endl;
  }
  

  *(G4CMPProcessUtils*)partitioner = *(G4CMPProcessUtils*)this;
  partitioner->SetVerboseLevel(verboseLevel);
  partitioner->UseVolume(aTrack.GetVolume());

  // Each charge carrier is independent, so it only gives back 0.5 times
  // the band gap. When the charge recombines, it may be in a different
  // location than its partner hole (due to a bias voltage), so combining
  // the two tracks is neither reasonable, nor The Geant4 Way.
  G4double Erecomb = (0.5*theLattice->GetBandGapEnergy()
		      + aTrack.GetKineticEnergy());

  partitioner->ComputePhononSampling(Erecomb);
  partitioner->DoPartition(0., Erecomb);
  partitioner->GetSecondaries(&aParticleChange);

  if (aParticleChange.GetNumberOfSecondaries() == 0) {	// Record energy release
    aParticleChange.ProposeNonIonizingEnergyDeposit(Erecomb);
  }

  aParticleChange.ProposeEnergy(0.);
  aParticleChange.ProposeTrackStatus(fStopAndKill);

  ClearNumberOfInteractionLengthLeft();		// All processes should do this!
  return &aParticleChange;
}


// Tracks should recombine if stopped, or if permanently below NTL threshold

G4bool G4CMPDriftRecombinationProcess::
ReadyToRecombine(const G4Track& aTrack) const {
  // SPECIAL: We use fStopAndKill to destroy tracks without recombining
  if (aTrack.GetTrackStatus() == fStopAndKill) return false;
    
  if (aTrack.GetTrackStatus() == fStopButAlive) {
    if (verboseLevel>1) G4cout << " track stopped." << G4endl;
    return true;
  }

  if (aTrack.GetStepLength() <= 0.) return false;	// Avoid reflections

  if (verboseLevel>1)
    G4cout << GetProcessName() << "::ReadyToRecombine?" << G4endl;

  // Non-physical recombination scale means no bulk recombination
  G4double mfpScale = G4CMPConfigManager::GetRecombinationScale();
  if (mfpScale <= 0. || mfpScale == DBL_MAX) {
    if (verboseLevel>1) G4cout << " bulk recombination suppressed." << G4endl;
    return false;
  }

  // Recombine now if no NTL emission expected before hitting surface
  G4bool noNTL = !LukeBeforeSurface(aTrack);
  if (verboseLevel>1) {
    G4cout << " " << (noNTL?"No ":"") << "NTL emission before surface."
	   << G4endl;
  }

  // Recombine now if no energy gain expected before hitting surface
  G4bool noGain = !EnergyGainToSurface(aTrack);
  if (verboseLevel>1) {
    G4cout << " " << (noGain?"No ":"") << "energy gained before surface."
	   << G4endl;
  }

  return (noNTL && noGain);
}


// Distance and direction to surface, along track or along field
// NOTE: Assumes "turn around" if electric field can redirect track

G4ThreeVector G4CMPDriftRecombinationProcess::
VectorToSurface(const G4Track& aTrack) const {
  // Use track touchable to create wrapper for G4VSolid interface
  G4CMPSolidUtils sutil(aTrack.GetTouchable(),verboseLevel,"Recombination");

  G4ThreeVector vSurf = GetGlobalMomentum(aTrack).unit();

  G4ThreeVector Efield = G4CMP::GetFieldAtPosition(aTrack);
  Efield *= aTrack.GetDynamicParticle()->GetCharge();

  if (Efield.mag() > 0.) vSurf = GetAcceleration(Efield).unit();

  G4double vDist = sutil.GetDistanceToSolid(aTrack.GetPosition(), vSurf);

  if (verboseLevel>2) {
    G4cout << " VectorToSurface: vDist " << vDist/mm << " mm"
	   << "   along " << vSurf << G4endl;
  }

  return vDist*vSurf;
}


// Compute maximum energy gain along current trajectory toward surface
// NOTE: Will include "turn around" if electric field is in use

G4bool G4CMPDriftRecombinationProcess::
EnergyGainToSurface(const G4Track& aTrack) const {
  G4ThreeVector Efield = G4CMP::GetFieldAtPosition(aTrack);
  if (Efield.mag() <= 0.) return false;

  // Energy gained from current position due to voltage bias
  G4ThreeVector vSurf = VectorToSurface(aTrack);
  G4double Egain = fabs(vSurf*Efield);

  // Minimum energy for NTL emission is kinetic energy at Vsound
  G4double Eluke = GetRateModel()->Threshold(1e-9*eV);
  G4double Etrack = GetKineticEnergy(aTrack);

  // Could track exceed NTL threshold before reaching surface?
  G4bool hasGain = (Etrack+Egain > Eluke);

  if (verboseLevel>2) {
    G4cout << " To surface " << vSurf.mag()/mm << " mm"
	   << " Egain " << Egain/eV << " eV: Egain+Etrack "
	   << (hasGain ? "exceeds" : "does not exceeed")
	   << " Eluke " << Eluke/eV << " eV" << G4endl;
  }

  return hasGain;
}


// Compare MFP for next NTL emission with distance to surface

G4bool G4CMPDriftRecombinationProcess::
LukeBeforeSurface(const G4Track& aTrack) const {
  G4double lukeMFP = GetMFPfromRate(aTrack);	  // Registered NTL RateModel
  G4ThreeVector vSurf = VectorToSurface(aTrack);  // Distance to track end

  // Rough guess at when next NTL emission might occur, with scale factor
  G4double mfpScale = G4CMPConfigManager::GetRecombinationScale();
  G4bool couldNTL = (mfpScale*lukeMFP < vSurf.mag());

  if (verboseLevel>2) {
    G4cout << " LukeBeforeSurface dist " << vSurf.mag()/mm << " mm"
	   << " NTL MFP " << lukeMFP/mm << " mm" << G4endl
	   << "   " << (couldNTL?"could":"probably won't")
	   << " emit NTL phonon" << G4endl;
  } 

  return couldNTL;
}


// Compute acceleration direction from Efield
// Apply H-V transform here for electrons

G4ThreeVector G4CMPDriftRecombinationProcess::
GetAcceleration(const G4ThreeVector& Efield) const {
  if (verboseLevel>2)
    G4cout << GetProcessName() << "::GetAcceleration" << G4endl;

  G4ThreeVector accel = Efield;

  if (IsHole()) return accel;			// Holes are simple particles

  // Rotate force into and out of valley frame, applying Herring-Vogt transform
  G4int valleyIndex = GetCurrentValley();
  const G4RotationMatrix& nToV = theLattice->GetValley(valleyIndex);
  const G4RotationMatrix& vToN = theLattice->GetValleyInv(valleyIndex);

  RotateToLocalDirection(accel);
  theLattice->RotateToLattice(accel);
  accel.transform(nToV);			// Rotate to valley frame
  accel *= theLattice->GetMInvTensor();
  accel *= theLattice->GetElectronMass();
  accel.transform(vToN);			// Back to lattice
  theLattice->RotateToSolid(accel);		// Back to crystal frame
  RotateToGlobalDirection(accel);

  if (verboseLevel>2) {
    G4cout << " Efield along " << Efield.unit() << " has acceleration"
	   << " along " << accel.unit() << G4endl;
  }

  return accel;
}

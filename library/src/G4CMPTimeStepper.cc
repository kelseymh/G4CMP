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
// 20170905  Cache Luke and IV rate models in local LoadDataFromTrack()
// 20170908  Remove "/10." rescaling of field when computing steps
// 20170919  Use rate threshold interface to define alternate step lengths
// 20180831  Fix compiler warning with PostStepDoIt() arguments
// 20190906  Provide functions to externally set rate models, move process
//		lookup functionality to G4CMP(Track)Utils.
// 20200331  C. Stanford (G4CMP-195): Added charge trapping
// 20200331  G4CMP-196: Added impact ionization mean free path
// 20200426  G4CMP-196: Use static function in TrapIonization for MFPs
// 20200504  M. Kelsey (G4CMP-195):  Get trapping MFPs from process
// 20200520  "First report" flag must be thread-local.
// 20200804  Move field access to G4CMPFieldUtils

#include "G4CMPTimeStepper.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4CMPDriftTrackInfo.hh"
#include "G4CMPFieldUtils.hh"
#include "G4CMPDriftTrappingProcess.hh"
#include "G4CMPDriftTrapIonization.hh"
#include "G4CMPGeometryUtils.hh"
#include "G4CMPTrackUtils.hh"
#include "G4CMPUtils.hh"
#include "G4CMPVProcess.hh"
#include "G4CMPVScatteringRate.hh"
#include "G4DynamicParticle.hh"
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
  : G4CMPVDriftProcess("G4CMPTimeStepper", fTimeStepper), lukeRate(nullptr),
    ivRate(nullptr), trappingLength(0.), trapIonLength(0.) {;}

G4CMPTimeStepper::~G4CMPTimeStepper() {;}


// Get scattering rates from current track's processes

void G4CMPTimeStepper::LoadDataForTrack(const G4Track* aTrack) {
  G4CMPProcessUtils::LoadDataForTrack(aTrack);	// Common configuration

  // Get rate model for Luke phonon emission from process
  const G4CMPVProcess* lukeProc =
    dynamic_cast<G4CMPVProcess*>(G4CMP::FindProcess(aTrack,
						    "G4CMPLukeScattering"));
  lukeRate = lukeProc ? lukeProc->GetRateModel() : nullptr;

  // get rate model for intervalley scattering from process
  const G4CMPVProcess* ivProc =
    dynamic_cast<G4CMPVProcess*>(G4CMP::FindProcess(aTrack,
					    "G4CMPInterValleyScattering"));
  ivRate = ivProc ? ivProc->GetRateModel() : nullptr;

  // get charge trapping mean free path
  trappingLength =
    G4CMPDriftTrappingProcess::GetMeanFreePath(GetCurrentParticle());

  // get charge-trap ionization mean free path
  G4double eTrapIonMFP = G4CMPDriftTrapIonization::
    GetMeanFreePath(GetCurrentParticle(), G4CMPDriftElectron::Definition());
  G4double hTrapIonMFP = G4CMPDriftTrapIonization::
    GetMeanFreePath(GetCurrentParticle(), G4CMPDriftHole::Definition());
  trapIonLength = std::min(eTrapIonMFP, hTrapIonMFP);

  if (verboseLevel>1) {
    G4cout << "TimeStepper Found" 
	   << (lukeRate?" lukeRate":"")
	   << (ivRate?" ivRate":"") 
	   << (trappingLength?" trappingLength":"") 
	   << (trapIonLength?" trapIonLength":"") 
	   << G4endl;
  }

  // Adjust rate models to keep highest verbosity level
  if (lukeRate && lukeRate->GetVerboseLevel() < verboseLevel)
    const_cast<G4CMPVScatteringRate*>(lukeRate)->SetVerboseLevel(verboseLevel);

  if (ivRate && ivRate->GetVerboseLevel() < verboseLevel)
    const_cast<G4CMPVScatteringRate*>(ivRate)->SetVerboseLevel(verboseLevel);
}


// Compute fixed "minimum distance" to avoid accelerating past Luke or IV

G4double G4CMPTimeStepper::GetMeanFreePath(const G4Track& aTrack, G4double,
					   G4ForceCondition* cond) {
  if (verboseLevel == -1) ReportRates(aTrack);	// SPECIAL FLAG TO REPORT

  *cond = NotForced;

  // SPECIAL:  If no electric field, no need to limit steps
  if (G4CMP::GetFieldAtPosition(aTrack).mag() <= 0.) return DBL_MAX;

  // Evaluate different step lengths to avoid overrunning process thresholds
  G4double vtrk = GetVelocity(aTrack);
  G4double ekin = GetKineticEnergy(aTrack);

  // Get step length due to fastest process
  G4double rate0 = MaxRate(aTrack);
  G4double mfp0 = rate0>0. ? vtrk/rate0 : DBL_MAX;

  if (verboseLevel>1) {
    G4cout << "TS Vtrk " << vtrk/(m/s) << " m/s mfp0 " << mfp0/m << " m"
	   << G4endl;
  }

  // Find distance to Luke threshold
  G4double mfp1 = lukeRate ? EnergyStep(lukeRate->Threshold(ekin)) : DBL_MAX;
  if (mfp1 <= 1e-9*m) mfp1 = DBL_MAX;	// Avoid steps getting "too short"

  if (verboseLevel>1)
    G4cout << "TS Luke threshold mfp1 " << mfp1/m << " m" << G4endl;

  // Find distance to IV scattering threshold 
  G4double mfp2 = ivRate ? EnergyStep(ivRate->Threshold(ekin)) : DBL_MAX;
  if (mfp2 <= 1e-9*m) mfp2 = DBL_MAX;	// Avoid steps getting "too short"

  if (verboseLevel>1 && ivRate)
    G4cout << "TS IV threshold mfp2 " << mfp2/m << " m" << G4endl;

  // Find distance to impact ionization
  G4double mfp3 = trapIonLength;

  if (verboseLevel>1)
    G4cout << "TS trap ionization mfp3 " << mfp3/m << " m" << G4endl;

  // Find MFP for charge trapping
  G4double mfp4 = trappingLength;
  if (mfp4 <= 1e-9*m) mfp4 = DBL_MAX;	// Avoid steps getting "too short"

  if (verboseLevel>1 && trappingLength)
    G4cout << "TS trapping MFP mfp4 " << mfp4/m << " m" << G4endl;

  // Take shortest distance from above options
  G4double mfp = std::min({mfp0, mfp1, mfp2, mfp3, mfp4});

  if (verboseLevel) {
    G4cout << GetProcessName() << (IsElectron()?" elec":" hole")
	   << " MFP = " << mfp/m << " m" << G4endl;
  }

  return mfp;
}


// At end of step, recompute kinematics; important for electrons

G4VParticleChange* G4CMPTimeStepper::PostStepDoIt(const G4Track& aTrack,
						  const G4Step& /*aStep*/) {
  aParticleChange.Initialize(aTrack);

  // Adjust mass and kinetic energy using end-of-step momentum
  G4ThreeVector pfinal = GetGlobalMomentum(aTrack);
  FillParticleChange(GetValleyIndex(aTrack), pfinal);

  ClearNumberOfInteractionLengthLeft();		// All processes must do this!
  return &aParticleChange;
}


// Get maximum rate for other processes at given kinematics

G4double G4CMPTimeStepper::MaxRate(const G4Track& aTrack) const {
  G4double lrate = lukeRate ? lukeRate->Rate(aTrack) : 0.;
  G4double irate = ivRate ? ivRate->Rate(aTrack) : 0.;

  if (verboseLevel>2) {
    G4cout << "G4CMPTimeStepper::MaxRate luke " << lrate/hertz << " iv "
	   << irate/hertz << " Hz" << G4endl;
  }

  return std::max(lrate,irate);
}


// Get step length in E-field needed to reach specified energy

G4double G4CMPTimeStepper::EnergyStep(G4double Efinal) const {
  const G4Track* trk = GetCurrentTrack();

  G4double Emag = G4CMP::GetFieldAtPosition(*trk).mag();
  if (Emag <= 0.) return DBL_MAX;		// No field, no acceleration

  G4double Ekin = GetKineticEnergy(trk);
  if (Ekin > Efinal) return DBL_MAX;		// Already over threshold

  if (verboseLevel>1) {
    G4cout << "G4CMPTimeStepper::EnergyStep from " << Ekin/eV
	   << " to " << Efinal/eV << " eV" << G4endl;
  }

  // Add 20% rescaling to account for electron valley systematics
  return 1.2*(Efinal-Ekin)/Emag;
}


// Report Luke and IV rates for diagnostics

void G4CMPTimeStepper::ReportRates(const G4Track& aTrack) {
  static G4ThreadLocal G4bool first = true;
  if (first) {
    G4cout << "TSreport Process Chg E[eV] v[m/s] k[1/m] Rate[Hz]" << G4endl;
    first = false;
  }
  
  G4cout << "TSreport Luke"
	 << " " << aTrack.GetDynamicParticle()->GetCharge()/eplus
	 << " " << GetKineticEnergy(aTrack)/eV
	 << " " << GetVelocity(aTrack)/(m/s)
	 << " " << GetLocalWaveVector(aTrack).mag()*m
	 << " " << lukeRate->Rate(aTrack)/hertz << G4endl;
  
  if (IsElectron()) {
    G4cout << "TSreport IV"
	   << " " << aTrack.GetDynamicParticle()->GetCharge()/eplus
	   << " " << GetKineticEnergy(aTrack)/eV
	   << " " << GetVelocity(aTrack)/(m/s)
	   << " " << GetLocalWaveVector(aTrack).mag()*m
	   << " " << ivRate->Rate(aTrack)/hertz << G4endl;
  }      
}

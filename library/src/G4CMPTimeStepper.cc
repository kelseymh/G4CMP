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
// 20210923  Ensure that rate calculations are initialized for track
// 20220504  E. Michaud -- Change DBL_MAX to minStep, remove negative MFPs
// 20220729  M. Kelsey -- EnergyStep() has incorrect units for output: should
//		be delta(E)/(q*V).
// 20220730  Drop trapping processes, as they have built-in MFPs, and don't
//		need TimeStepper for energy-dependent calculation.

#include "G4CMPTimeStepper.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4CMPDriftTrackInfo.hh"
#include "G4CMPFieldUtils.hh"
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
    ivRate(nullptr) {;}

G4CMPTimeStepper::~G4CMPTimeStepper() {;}


// Get scattering rates from current track's processes

void G4CMPTimeStepper::LoadDataForTrack(const G4Track* aTrack) {
  G4CMPProcessUtils::LoadDataForTrack(aTrack);	// Common configuration

  // Get rate model for Luke phonon emission from process
  const G4CMPVProcess* lukeProc =
    dynamic_cast<G4CMPVProcess*>(G4CMP::FindProcess(aTrack,
						    "G4CMPLukeScattering"));
  lukeRate = lukeProc ? lukeProc->GetRateModel() : nullptr;
  if (lukeRate)
    const_cast<G4CMPVScatteringRate*>(lukeRate)->LoadDataForTrack(aTrack);

  // get rate model for intervalley scattering from process
  const G4CMPVProcess* ivProc =
    dynamic_cast<G4CMPVProcess*>(G4CMP::FindProcess(aTrack,
					    "G4CMPInterValleyScattering"));
  ivRate = ivProc ? ivProc->GetRateModel() : nullptr;
  if (ivRate) 
    const_cast<G4CMPVScatteringRate*>(ivRate)->LoadDataForTrack(aTrack);

  if (verboseLevel>1) {
    G4cout << "TimeStepper Found" 
	   << (lukeRate?" lukeRate":"")
	   << (ivRate?" ivRate":"") 
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
  G4double rate = MaxRate(aTrack);
  G4double mfpFast = rate>0. ? vtrk/rate : DBL_MAX;

  if (verboseLevel>1) {
    G4cout << "TS Vtrk " << vtrk/(m/s) << " m/s mfp0 " << mfpFast/m << " m"
	   << G4endl;
  }


  // Find distance to Luke threshold
  G4double mfpLuke = DistanceToThreshold(lukeRate, ekin);
  if (verboseLevel>1)
    G4cout << "TS Luke threshold mfpLuke " << mfpLuke/m << " m" << G4endl;

  // Find distance to IV scattering threshold 
  G4double mfpIV = DistanceToThreshold(ivRate, ekin);
  if (verboseLevel>1)
    G4cout << "TS IV threshold mfpIV " << mfpIV/m << " m" << G4endl;

  // Take shortest distance from above options
  G4double mfp = std::min({mfpFast, mfpLuke, mfpIV});

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


// Get distance due to energy gain to threshold (step) in rate

G4double G4CMPTimeStepper::DistanceToThreshold(const G4CMPVScatteringRate* rate,
					       G4double Estart) const {
  if (!rate) return DBL_MAX;		// Skip if no rate model

  // Avoid taking "too short" steps, which causes "stuck tracks"
  G4double MINstep = G4CMPConfigManager::GetMinStepScale();
  MINstep *= (IsElectron() ? theLattice->GetElectronScatter()
		: theLattice->GetHoleScatter());

  if (MINstep<0) MINstep = 1e-6*m;

  return std::max(EnergyStep(Estart, rate->Threshold(Estart)), MINstep);
}

// Get step length in E-field needed to reach specified energy

G4double G4CMPTimeStepper::EnergyStep(G4double Estart, G4double Efinal) const {
  if (Estart > Efinal) return DBL_MAX;		// Already over threshold

  const G4Track* trk = GetCurrentTrack();

  G4double EMmag = G4CMP::GetFieldAtPosition(*trk).mag();
  if (EMmag <= 0.) return DBL_MAX;		// No field, no acceleration

  if (verboseLevel>1) {
    G4cout << "G4CMPTimeStepper::EnergyStep from " << Estart/eV
	   << " to " << Efinal/eV << " eV : " 
	   << (Efinal-Estart)/(eplus*EMmag) << " mm" << G4endl;
  }

  // Add 20% rescaling to account for electron valley systematics
  return 1.2*(Efinal-Estart)/(eplus*EMmag);
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

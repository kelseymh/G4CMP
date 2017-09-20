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

#include "G4CMPTimeStepper.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4CMPDriftTrackInfo.hh"
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
  : G4CMPVDriftProcess("G4CMPTimeStepper", fTimeStepper), maxStep(0.),
    tempTrack(nullptr), lukeRate(nullptr), ivRate(nullptr) {;}

G4CMPTimeStepper::~G4CMPTimeStepper() {
  delete tempTrack;
}


// Get scattering rates from current track's processes

void G4CMPTimeStepper::LoadDataForTrack(const G4Track* aTrack) {
  G4CMPProcessUtils::LoadDataForTrack(aTrack);	// Common configuration

  // Use maximum step length scaled to "scattering length"
  // NOTE:  Use of MinStepScale double-counts with G4CMPVDriftProcess!
  maxStep = 0.05; /* G4CMPConfigManager::GetMinStepScale(); */
  maxStep *= (IsElectron() ? theLattice->GetElectronScatter()
	      : theLattice->GetHoleScatter());

  if (verboseLevel>2)
    G4cout << "TimeStepper maxStep " << maxStep/mm << " mm for "
	   << aTrack->GetDefinition()->GetParticleName() << G4endl;

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

  // Take shortest distance or minimum step length
  G4double mfp = std::min(std::min(mfp0, mfp1), mfp2);
  if (maxStep > 0.) mfp = std::min(mfp, maxStep);

  if (verboseLevel) {
    G4cout << GetProcessName() << (IsElectron()?" elec":" hole")
	   << " MFP = " << mfp/m << " m" << G4endl;
  }

  return mfp;
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


// Duplicate track information for step evaluation

void G4CMPTimeStepper::CopyTrack(const G4Track& aTrack) {
  delete tempTrack;
  tempTrack = new G4Track(aTrack);
  tempTrack->SetTouchableHandle(aTrack.GetTouchableHandle());
  tempTrack->SetNextTouchableHandle(aTrack.GetNextTouchableHandle());
  tempTrack->SetOriginTouchableHandle(aTrack.GetOriginTouchableHandle());

  G4int infoID = G4CMPConfigManager::GetPhysicsModelID();
  auto info = G4CMP::GetTrackInfo<G4CMPDriftTrackInfo>(aTrack);

  G4CMPDriftTrackInfo* tempInfo =
    new G4CMPDriftTrackInfo(info->Lattice(), info->ValleyIndex());

  tempTrack->SetAuxiliaryTrackInformation(infoID, tempInfo);
}


// Report Luke and IV rates for diagnostics

void G4CMPTimeStepper::ReportRates(const G4Track& aTrack) {
  static G4bool first = true;
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

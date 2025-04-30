/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPVProcess.cch
/// \brief Top-level base class for all G4CMP physics processes.  Only
///        discrete (post-step) processes are supported in G4CMP.  This
///	   base class provides access to all of the "ProcessUtils"
///	   functionaltiy via multiple inheritance.  Concrete processes
///	   don't need to do anything special.
//
// $Id$
//
// 20170601  New abstract base class for all G4CMP processes
// 20170802  Add registration of external scattering rate (MFP) model
// 20170815  Call through to scattering-rate LoadDataForTrack()
// 20190906  Bug fix in UseRateModel(), check for good pointer, not null;
//		Add function to initialize rate model after LoadDataForTrack
// 20210915  Change diagnostic output to verbose=3 or higher.
// 20250430  Add ability to register different physics models, which will
//	       be called by (new) base implementation of PostStepDoIt().

#include "G4CMPVProcess.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPVPhysicsModel.hh"
#include "G4CMPVScatteringRate.hh"
#include "G4ForceCondition.hh"
#include "G4SystemOfUnits.hh"


// Constructor and destructor

G4CMPVProcess::G4CMPVProcess(const G4String& processName,
			     G4CMPProcessSubType stype)
  : G4VDiscreteProcess(processName, fPhonon), G4CMPProcessUtils(),
    rateModel(0), physicsModel(0) {
  verboseLevel = G4CMPConfigManager::GetVerboseLevel();
  SetProcessSubType(stype);
}

G4CMPVProcess::~G4CMPVProcess() {
  delete rateModel; rateModel=0;
  delete physicsModel; physicsModel=0;
}


// Register utility class for computing scattering rate for MFP
// NOTE:  Takes ownership of model for deletion; deletes any previous version

void G4CMPVProcess::UseRateModel(G4CMPVScatteringRate* model) {
  if (model == rateModel) return;		// Nothing to change

  if (rateModel) delete rateModel;		// Avoid memory leaks!
  rateModel = model;

  // Ensure that rate model is syncronized with process state
  ConfigureRateModel();
}


// Register utility class for implementing particular scattering model
// NOTE:  Takes ownership of model for deletion; deletes any previous version

void G4CMPVProcess::UsePhysicsModel(G4CMPVPhysicsModel* model) {
  if (model == physicsModel) return;		// Nothing to change

  if (physicsModel) delete physicsModel;	// Avoid memory leaks!
  physicsModel = model;

  model->SetProcess(this);			// Ensure updated configuration
  model->SetVerboseLevel(verboseLevel);
}


// Configuration

void G4CMPVProcess::StartTracking(G4Track* track) {
  G4VProcess::StartTracking(track);	// Apply base class actions
  LoadDataForTrack(track);
  ConfigureRateModel();
}

void G4CMPVProcess::EndTracking() {
  G4VProcess::EndTracking();		// Apply base class actions
  ReleaseTrack();
  if (rateModel) rateModel->ReleaseTrack();
}

void G4CMPVProcess::ConfigureRateModel() {
  if (!rateModel) return;

  rateModel->SetVerboseLevel(verboseLevel);
  if (GetCurrentTrack()) rateModel->LoadDataForTrack(GetCurrentTrack());
}


// Use physics model to generate final state kinematics
// Concrete processes which don't use a model will override this

G4VParticleChange* 
G4CMPVProcess::PostStepDoIt(const G4Track& track, const G4Step& step) {
  return (physicsModel ? physicsModel->PostStepDoIt(track,step)
	  : G4VDiscreteProcess::PostStepDoIt(track,step));
}

// Compute MFP using track velocity and scattering rate

G4double G4CMPVProcess::GetMeanFreePath(const G4Track& aTrack, G4double,
					G4ForceCondition* condition) {
  *condition = (rateModel && rateModel->IsForced()) ? Forced : NotForced;

  G4double rate = rateModel ? rateModel->Rate(aTrack) : 0.;
  G4double vtrk = IsChargeCarrier() ? GetVelocity(aTrack) : aTrack.GetVelocity();
  G4double mfp  = rate>0. ? vtrk/rate : DBL_MAX;

  if (verboseLevel>2) {
    G4cout << GetProcessName() << " rate = " << rate/hertz << " Hz"
	   << " Vtrk = " << vtrk/(m/s) << " m/s"
	   << " MFP = " << mfp/m << " m" << G4endl;
  }

  return mfp;
}

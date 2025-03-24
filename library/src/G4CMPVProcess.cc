/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPVProcess.cch
/// \brief Implementation of the G4CMPVProcess base class
//
// $Id$
//
// 20170601  New abstract base class for all G4CMP processes
// 20170802  Add registration of external scattering rate (MFP) model
// 20170815  Call through to scattering-rate LoadDataForTrack()
// 20190906  Bug fix in UseRateModel(), check for good pointer, not null;
//		Add function to initialize rate model after LoadDataForTrack
// 20210915  Change diagnostic output to verbose=3 or higher.

#include "G4CMPVProcess.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPVScatteringRate.hh"
#include "G4ForceCondition.hh"
#include "G4SystemOfUnits.hh"
#include "G4CMPTrackUtils.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"

// Constructor and destructor

G4CMPVProcess::G4CMPVProcess(const G4String& processName,
			     G4CMPProcessSubType stype)
  : G4VDiscreteProcess(processName, fPhonon), G4CMPProcessUtils(), G4CMPSCUtils(),
    rateModel(0) {
  verboseLevel = G4CMPConfigManager::GetVerboseLevel();
  SetProcessSubType(stype);
}

G4CMPVProcess::~G4CMPVProcess() {
  delete rateModel; rateModel=0;
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


//This logic block needs to be run in every process's GetMeanFreePath in order to let it know that the lattice has changed.
//Otherwise the different processes (which inherit from individual/separate instances of the G4CMPProcessUtils class) will
//see different lattices.
G4bool G4CMPVProcess::UpdateMeanFreePathForLatticeChangeover(const G4Track& aTrack)
{
  if( verboseLevel > 5 ){
    G4cout << "---------- G4CMPVProcess::UpdateMeanFreePathForLatticeChangeover ----------" << G4endl;
    G4cout << "UMFPFLC Function Point A | loading data for track after lattice changeover, process: " << this->GetProcessName() << G4endl;
    G4cout << "UMFPFLC Function Point A | Here, track length: " << aTrack.GetTrackLength() << G4endl;
    G4cout << "UMFPFLC Function Point A | Current lattice a la lattice manager: " << G4LatticeManager::GetLatticeManager()->GetLattice(aTrack.GetVolume()) << ", while this->theLattice: " << this->theLattice << G4endl;
  }
    
  //Always do a check to see if the current lattice stored in this process is equal to the one that represents
  //the volume that we're in. Note that we can't do this with the "GetLattice()" and "GetNextLattice()" calls
  //here because at this point in the step, the pre- and post-step points both point to the same volume. Since
  //GetMeanFreePath is run at the beginning, I think the point at which a boundary interaction is assessed comes
  //later (hence why we can use that info in PostStepDoIts but not here.) Adding a statement about track length here,
  //since it seems that when a particle spawns it doesn't necessarily trigger this block, and I think we want it to.
  if( (((this->theLattice) && G4LatticeManager::GetLatticeManager()->GetLattice(aTrack.GetVolume())) &&
       (this->theLattice != G4LatticeManager::GetLatticeManager()->GetLattice(aTrack.GetVolume()))) ||
      aTrack.GetTrackLength() == 0.0 ){

    //Debugging
    if( verboseLevel > 5 ){
      G4cout << "UMFPFLC Function Point B | the step length associated with this is " << aTrack.GetStep()->GetStepLength() << G4endl;
    }
    
    
    //REL noting that if physical lattices are not 1:1 with volumes, something may get broken here... Should check a scenario of segmented SC...
    
    this->LoadDataForTrack(&aTrack);
    if(rateModel) rateModel->LoadDataForTrack(&aTrack);

    //Debugging
    if( verboseLevel > 5 ){
      G4cout << "UMFPFLC Function Point C | Successfully changed over to a new lattice for process " << this->GetProcessName() << G4endl;
    }
    return true;
    
  }
  //Debugging
  if( verboseLevel > 5 ){    
    G4cout << "UMFPFLC Function Point D | Did not successfully change over to a new lattice for process " << this->GetProcessName() << G4endl;
  }
  return false;
}

//This is meant to update superconductor info for both the process and the rate info if we move into a new lattice
void G4CMPVProcess::UpdateSCAfterLatticeChange()
{
  //Debugging
  if( verboseLevel > 5 ){
    G4cout << "---------- G4CMPVProcess::UpdateSCAfterLatticeChange ----------" << G4endl;
  }
 
  //First, determine if the new lattice is a SC. If not, then set the SCUtils info to null for this process  
  if( (this->theLattice)->GetSCDelta0() <= 0 ){
    this->SetCurrentSCInfoToNull();
    if(rateModel) rateModel->SetCurrentSCInfoToNull();
    return;
  }

  //If it is a SC, then we should update the SC information for the SC utils class within the base of this
  //and the base of the rate model. Also, handle the checking/updating of the lookup tables to be used for each
  //SC.
  this->LoadLatticeInfoIntoSCUtils(this->theLattice);
  if( rateModel ){
    rateModel->LoadLatticeInfoIntoSCUtils(this->theLattice); 
    rateModel->UpdateLookupTable(this->theLattice);
  }
}



// Compute MFP using track velocity and scattering rate

G4double G4CMPVProcess::GetMeanFreePath(const G4Track& aTrack, G4double,
					G4ForceCondition* condition) {
  //Debugging
  if( verboseLevel > 5 ){
    G4cout << "---------- G4CMPVProcess::GetMeanFreePath ----------" << G4endl;
  }

  //Update lattice information within the process utils 
  if(UpdateMeanFreePathForLatticeChangeover(aTrack)){
    UpdateSCAfterLatticeChange();
  }

  
  *condition = (rateModel && rateModel->IsForced()) ? Forced : NotForced;


  G4double rate = rateModel ? rateModel->Rate(aTrack) : 0.;
  G4double vtrk = IsChargeCarrier() ? GetVelocity(aTrack) : aTrack.GetVelocity();
  G4double mfp  = rate>0. ? vtrk/rate : DBL_MAX;

  //Very important debugging
  if( verboseLevel > 5 ){    
    G4cout << "GMFP Function Point A | In getMFP, rate of " << this->GetProcessName() << " is: " << rate << G4endl;
  }
  if (verboseLevel>2) {
    G4cout << GetProcessName() << " rate = " << rate/hertz << " Hz"
	   << " Vtrk = " << vtrk/(m/s) << " m/s"
	   << " MFP = " << mfp/m << " m" << G4endl;
  }

  return mfp;
}

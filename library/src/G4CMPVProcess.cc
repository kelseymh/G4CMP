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
// 20250905  G4CMP-500 -- Now using a fundamental SC parameter (i.e. not
//              the gap0) to determine if we're in a superconducting volume

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
  : G4VDiscreteProcess(processName, fPhonon), G4CMPProcessUtils(),
    G4CMPSCUtils(),
    rateModel(0) {
  verboseLevel = G4CMPConfigManager::GetVerboseLevel(); 
  //G4cout << "In G4CMPVProcess, setting verbose level to: " << verboseLevel
  //<< G4endl;
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


//This logic block needs to be run in every process's GetMeanFreePath in order
//to let it know that the lattice has changed. Otherwise the different
//processes (which inherit from individual/separate instances of the
//G4CMPProcessUtils class) will see different lattices.
G4bool G4CMPVProcess::
UpdateMeanFreePathForLatticeChangeover(const G4Track& aTrack) {
  if (verboseLevel > 5) {
    G4cout << "-- G4CMPVProcess::UpdateMeanFreePathForLatticeChangeover --"
           << G4endl;
    G4cout << "UMFPFLC Function Point A | loading data for track after lattice "
           << "changeover, process: " << this->GetProcessName() << G4endl;
    G4cout << "UMFPFLC Function Point A | Here, track length: "
           << aTrack.GetTrackLength() << G4endl;
    G4cout << "UMFPFLC Function Point A | Current lattice a la lattice "
           << "manager: "
           << G4LatticeManager::GetLatticeManager()->GetLattice(aTrack.GetVolume())
           << ", while this->theLattice: " << this->theLattice << G4endl;
  }
    
  //Always do a check to see if the current lattice stored in this process is
  //equal to the one that represents the volume that we're in. Note that we
  //can't do this with the "GetLattice()" and a hypothetical "GetNextLattice()" 
  //because at this point in the step, the pre- and post-step points both
  //point to the same volume. Since GetMeanFreePath is run at the beginning, I
  //think the point at which a boundary interaction is assessed comes later
  //(hence why we can use that info in PostStepDoIts but not here.) Adding a
  //statement about track length here, since it seems that when a particle
  //spawns it doesn't necessarily trigger this block, and I think we want it to.
  if ((((this->theLattice) && G4LatticeManager::GetLatticeManager()->GetLattice(aTrack.GetVolume())) &&
       (this->theLattice != G4LatticeManager::GetLatticeManager()->GetLattice(aTrack.GetVolume()))) ||
      aTrack.GetTrackLength() == 0.0) {

    //Debugging
    if (verboseLevel > 5) {
      G4cout << "UMFPFLC Function Point B | the step length associated with "
             << "this is " << aTrack.GetStep()->GetStepLength() << G4endl;
    }
        
    //Noting here that since LoadDataForTrack updates the momentum based on
    //the current lattice (which in a reflection step is the "turnaround
    //lattice", we need to provide an opt-out for that momentum recalculation.
    //In some cases where the existing k-vector is shallow, that momentum
    //recalculation can take the existing k-vector and turn it around and
    //launch it into the new "turnaround" lattice in something that looks like
    //a transmission, not a reflection. So for phonons I just add in a
    //conditional to the LoadDataForTrack function that ensures it doesn't
    //attempt a recalculation of the momentum vector when loaddatafromtrack is
    //run from this specific location (i.e. the pre-step volume differs from
    //this->theLattice). We note that for actual, true transmissions, the
    //setting of the new k-vector and momentum in the new lattice should happen
    //*within* the doTransmission function at the end of the prior step. I
    //think that's the only case where this override could mess things up, and
    //it should be covered in that function already.
    this->LoadDataForTrack(&aTrack,true);
    if (rateModel) rateModel->LoadDataForTrack(&aTrack,true);

    //Debugging
    if (verboseLevel > 5) {
      G4cout << "UMFPFLC Function Point C | Successfully changed over to a new "
             << "lattice for process " << this->GetProcessName() << G4endl;
    }
    return true;
  }
  //Debugging
  if (verboseLevel > 5) {    
    G4cout << "UMFPFLC Function Point D | Did not successfully change over to "
           << "a new lattice for process " << this->GetProcessName() << G4endl;
  }
  return false;
}

//This is meant to update superconductor info for both the process and the rate
//info if we move into a new lattice
void G4CMPVProcess::UpdateSCAfterLatticeChange() {

  //Debugging
  if (verboseLevel > 5) {
    G4cout << "-- G4CMPVProcess::UpdateSCAfterLatticeChange --" << G4endl;
  }
 
  //First, determine if the new lattice is a SC. We do this via the QP relax-
  //ation time since the gap0 is not a fundamental/intrinsic property of the
  //material (since it's height-dependent).
  if ((this->theLattice)->GetSCTau0qp() == DBL_MAX) {
    this->SetCurrentSCInfoToNull();
    if (rateModel) rateModel->SetCurrentSCInfoToNull();
    return;
  }

  //If it is a SC, then we should update the SC information for the SC utils
  //class within the base of this and the base of the rate model. Also, handle
  //the checking/updating of the lookup tables to be used for each SC.
  this->LoadLatticeInfoIntoSCUtils(this->theLattice);
  if (rateModel) {
    rateModel->LoadLatticeInfoIntoSCUtils(this->theLattice); 
    rateModel->UpdateLookupTable(this->theLattice);
  }
}

// Compute MFP using track velocity and scattering rate
G4double G4CMPVProcess::GetMeanFreePath(const G4Track& aTrack, G4double,
					G4ForceCondition* condition) {    
  //Debugging
  if (verboseLevel > 5) {
    G4cout << "-- G4CMPVProcess::GetMeanFreePath --" << G4endl;
  }

  //Update lattice information within the process utils 
  if (UpdateMeanFreePathForLatticeChangeover(aTrack)) {
    UpdateSCAfterLatticeChange();
  }

  
  *condition = (rateModel && rateModel->IsForced()) ? Forced : NotForced;


  G4double rate = rateModel ? rateModel->Rate(aTrack) : 0.;
  G4double vtrk = IsChargeCarrier() ? GetVelocity(aTrack) : aTrack.GetVelocity();
  G4double mfp  = rate>0. ? vtrk/rate : DBL_MAX;

  //Very important debugging
  if (verboseLevel > 5) {    
    G4cout << "GMFP Function Point A | In getMFP, rate of "
	   << this->GetProcessName() << " is: " << rate << G4endl;
  }
  if (verboseLevel>2) {
    G4cout << GetProcessName() << " rate = " << rate/hertz << " Hz"
	   << " Vtrk = " << vtrk/(m/s) << " m/s"
	   << " MFP = " << mfp/m << " m" << G4endl;
  }
  return mfp;
}

/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPQPRecombinationProcess.cc
/// \brief Implementation of the G4CMPQPRecombinationProcess class
//
// $Id$
//

#include "G4CMPQPRecombinationProcess.hh"
#include "G4CMPQPRecombinationRate.hh"
#include "G4CMPSCUtils.hh"
#include "G4QP.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VParticleChange.hh"
#include "G4RandomDirection.hh"
#include "G4CMPUtils.hh"
#include "G4LatticePhysical.hh"
#include "G4CMPSecondaryUtils.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Constructor and destructor 
G4CMPQPRecombinationProcess::G4CMPQPRecombinationProcess(const G4String& aName)
  : G4VQPProcess(aName,fQPRecombinationProcess) {
  UseRateModel(new G4CMPQPRecombinationRate);
}

G4CMPQPRecombinationProcess::~G4CMPQPRecombinationProcess() {
}


void G4CMPQPRecombinationProcess::SetVerboseLevel(G4int vb) {
  verboseLevel = vb;
}

G4VParticleChange* G4CMPQPRecombinationProcess::
PostStepDoIt(const G4Track& aTrack,const G4Step& aStep) {

  //Debugging
  if (verboseLevel > 5) {
    G4cout << "-- G4CMPQPRecombinationProcess::PostStepDoIt --" << G4endl;
  }
  
  aParticleChange.Initialize(aTrack);
  
  //Pseudocode
  //1. Determine if we're on a boundary surface. If we are, kill the event --
  //   if the code is working properly, this should basically never happen...
  G4StepPoint * postStepPoint = aStep.GetPostStepPoint();
  if (postStepPoint->GetStepStatus() == fGeomBoundary ||
      postStepPoint->GetStepStatus() == fWorldBoundary) {
    G4ExceptionDescription msg;
    msg << "For some reason we're running post-step do it for the Bogoliubov "
	<< "Recombination process and we find ourselves on a boundary. Since "
	<< "the QPs are artificially told to move attometers by design, I'm "
	<< "not expecting this to happen.";
    G4Exception("G4CMPQPRecombinationProcess::PostStepDoIt",
		"QPRecombination001",EventMustBeAborted,msg);
    return &aParticleChange;		
  }

  //2. Determine if this QP is actually going to produce a phonon upon
  //   recombination -- to conserve energy, we have to restrict ourselves to
  //   only producing new 2Delta phonons from half of the QPs, since we're
  //   assuming recombination happens with a bath of ambient QPs that already
  //   exist independently of the scatter. Here, since QPs in any scenario
  //   other than a manual creation will be produced in pairs which should
  //   have adjacent track IDs (always), we only produce phonons if the track
  //   ID is even. So here, let's find our track id.
  G4int trackID = aTrack.GetTrackID();

  //2. If the track is odd, then create recombination phonon
  if (trackID % 2 == 1) {

    //Debugging
    if (verboseLevel > 5) {
      G4cout << "PSDI Function Point A | testing logic for recombination. "
	     << "TrackID, " << trackID << ", is odd, so a phonon is generated."
	     << G4endl;
    }
    
    //Make the energy of the new photon something simplistic for now: gap
    //energy plus energy of this QP.
    G4double newPhonEnergy = fGapEnergy + GetKineticEnergy(aTrack);
    GenerateRecombinationPhonon(newPhonEnergy,aTrack,aStep);    
  }
  
  //3. Otherwise, recombine and just kill it
  else {

    //Debugging
    if (verboseLevel > 5) {
      G4cout << "PSDI Function Point B | testing logic for recombination. "
	     << "TrackID, " << trackID << ", is even, and so no phonons are "
	     << "generated." << G4endl;
    }   
    aParticleChange.ProposeEnergy(0.);
    aParticleChange.ProposeTrackStatus(fStopAndKill);
  }
  
  //Don't have a "clear" process because the particle is dead at the end of
  //this function

  //6. Return the particle change
  return &aParticleChange;
}


G4bool
G4CMPQPRecombinationProcess::IsApplicable(const G4ParticleDefinition& aPD) {
  
  // Allow all phonon types, because type is changed during tracking
  return G4VQPProcess::IsApplicable(aPD);
}


// Take the phonon energy and produce a recombination phonon from this.
void G4CMPQPRecombinationProcess::
GenerateRecombinationPhonon(G4double phonEnergy,const G4Track& aTrack,
			    const G4Step& /*aStep*/) {

  //Debugging
  if (verboseLevel > 5) {
    G4cout << "-- G4CMPQPRecombinationProcess::GenerateRecombinationPhonon --"
	   << G4endl;
  }
  
  //Now create the phonon
  G4int mode = G4CMP::ChoosePhononPolarization(theLattice->GetLDOS(),
					       theLattice->GetSTDOS(),
					       theLattice->GetFTDOS());    
  G4ThreeVector dir1 = G4RandomDirection();    
  G4Track* sec1 = G4CMP::CreatePhonon(aTrack,mode,dir1,phonEnergy,
				      aTrack.GetGlobalTime(),
				      aTrack.GetPosition());

  //Debugging
  if (verboseLevel > 5) {
    G4cout << "GRP Function Point A | energy of recombination phonon: "
	   << phonEnergy << ", whereas twice the gap is: " << 2.0*fGapEnergy
	   << G4endl;
  }
  
  //Check to make sure the secondary was actually produced
  if (!sec1) {
    G4Exception("G4CMPQPRecombinationProcess::GenerateRecombinationPhonon",
		"QPRecombination002",
		JustWarning,
		"Error creating secondaries");
    return;
  }
  
  //Set the number of secondaries to 1 and add the secondary. Kill the QP.
  aParticleChange.SetNumberOfSecondaries(1);
  aParticleChange.AddSecondary(sec1);
  aParticleChange.ProposeEnergy(0.);
  aParticleChange.ProposeTrackStatus(fStopAndKill);
}
									

//Pass-through to G4CMPVProcess class
G4double G4CMPQPRecombinationProcess::GetMeanFreePath(const G4Track& trk,
						      G4double prevstep,
						      G4ForceCondition* cond) {
  //Debugging
  if (verboseLevel > 5) {
    G4cout << "-- G4CMPQPRecombinationProcess::GetMeanFreePath --" << G4endl;
  }
  
  G4double mfpBase = G4CMPVProcess::GetMeanFreePath(trk,prevstep,cond);

  //Debugging
  if (verboseLevel > 5) {
    G4cout << "GMFP Function Point A | mean free path in "
	   << "QPRecombinationProcess: " << mfpBase << ", with nMFPsLeft: "
	   << GetNumberOfInteractionLengthLeft() << G4endl;
  }

  return mfpBase;
}

/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPBogoliubovQPRecombinationProcess.cc
/// \brief Implementation of the G4CMPBogoliubovQPRecombinationProcess class
//
// $Id$
//

#include "G4CMPBogoliubovQPRecombinationProcess.hh"
#include "G4CMPBogoliubovQPRecombinationRate.hh"
#include "G4CMPSCUtils.hh"
#include "G4BogoliubovQP.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VParticleChange.hh"
#include "G4RandomDirection.hh"
#include "G4CMPUtils.hh"
#include "G4LatticePhysical.hh"
#include "G4CMPSecondaryUtils.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Constructor and destructor 
G4CMPBogoliubovQPRecombinationProcess::G4CMPBogoliubovQPRecombinationProcess(const G4String& aName)
  : G4VBogoliubovQPProcess(aName,fQPRecombinationProcess)
{  
  UseRateModel(new G4CMPBogoliubovQPRecombinationRate);
}

G4CMPBogoliubovQPRecombinationProcess::~G4CMPBogoliubovQPRecombinationProcess()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4CMPBogoliubovQPRecombinationProcess::SetVerboseLevel(G4int vb) {
  verboseLevel = vb;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4VParticleChange* G4CMPBogoliubovQPRecombinationProcess::PostStepDoIt(const G4Track& aTrack,
								       const G4Step& aStep)
{
  G4cout << "REL in BogoliubovQPRecombination process poststepdoit." << G4endl;
  
  aParticleChange.Initialize(aTrack);
  
  //Pseudocode
  //1. Determine if we're on a boundary surface. If we are, kill the event -- if the code is working properly, this should
  //   basically never happen...
  G4StepPoint * postStepPoint = aStep.GetPostStepPoint();
  if (postStepPoint->GetStepStatus() == fGeomBoundary ||
      postStepPoint->GetStepStatus() == fWorldBoundary) {
    G4ExceptionDescription msg;
    msg << "For some reason we're running post-step do it for the Bogoliubov Recombination process and we find ourselves on a boundary. Since the QPs are artificially told to move attometers by design, I'm not expecting this to happen.";
    G4Exception("G4CMPBogoliubovQPRecombinationProcess::PostStepDoIt", "BogoliubovQPRecombination001",EventMustBeAborted,msg);
    return &aParticleChange;		
  }

  //2. Determine if this QP is actually going to produce a phonon upon recombination -- to conserve energy,
  //   we have to restrict ourselves to only producing new 2Delta phonons from half of the QPs, since we're
  //   assuming recombination happens with a bath of ambient QPs that already exist independently of the scatter
  //   Here, we cheat and use velocity, which is set up to be either 1E-18 or 2E-18 -- velocity is useless anyway because
  //   our QP transport is going to have to be a special thing anyway.
  G4double vel = aTrack.GetVelocity() / (CLHEP::m / CLHEP::s);

  //2. If the velocity is 1E-18, then recombine and produce a phonon. Probably a better way to sense this
  if( vel / 1E-18 < 1.5 ){

    G4cout << "REL testing logic for recombination. Vel: " << vel << G4endl;
    
    //Make the energy of the new photon something simplistic for now: gap energy plus energy of this QP.
    G4double newPhonEnergy = fGapEnergy + GetKineticEnergy(aTrack);
    GenerateRecombinationPhonon(newPhonEnergy,aTrack,aStep);    
  }
  
  //3. Otherwise, recombine and just kill it
  else{
    aParticleChange.ProposeEnergy(0.);
    aParticleChange.ProposeTrackStatus(fStopAndKill);
  }

  //Don't have a "clear" process because the particle is dead at the end of this function
  
  //6. Return the particle change
  return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4bool G4CMPBogoliubovQPRecombinationProcess::IsApplicable(const G4ParticleDefinition& aPD) {
  // Allow all phonon types, because type is changed during tracking
  return G4VBogoliubovQPProcess::IsApplicable(aPD);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Take the phonon energy and produce a recombination phonon from this.
void G4CMPBogoliubovQPRecombinationProcess::GenerateRecombinationPhonon(G4double phonEnergy,
									const G4Track& aTrack,
									const G4Step& aStep)
{
  //Now create the phonon
  G4int mode = G4CMP::ChoosePhononPolarization(theLattice->GetLDOS(), theLattice->GetSTDOS(),theLattice->GetFTDOS());    
  G4ThreeVector dir1 = G4RandomDirection();    
  G4Track* sec1 = G4CMP::CreatePhonon(aTrack,mode,dir1,phonEnergy,aTrack.GetGlobalTime(),aTrack.GetPosition());

  G4cout << "REL energy of recombination phonon: " << phonEnergy << ", whereas twice the gap is: " << 2.0*fGapEnergy << G4endl;


  
  //Check to make sure the secondary was actually produced
  if (!sec1) {
    G4Exception("G4CMPBogoliubovQPRecombinationProcess::GenerateRecombinationPhonon", "BogoliubovQPRecombination002",
		JustWarning, "Error creating secondaries");
    return;
  }
  
  //Set the number of secondaries to 1 and add the secondary. Kill the QP.
  aParticleChange.SetNumberOfSecondaries(1);
  aParticleChange.AddSecondary(sec1);
  aParticleChange.ProposeEnergy(0.);
  aParticleChange.ProposeTrackStatus(fStopAndKill);
}
									
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//Pass-through to G4CMPVProcess class
G4double G4CMPBogoliubovQPRecombinationProcess::GetMeanFreePath(const G4Track& trk, G4double prevstep, G4ForceCondition* cond)
{
  return G4CMPVProcess::GetMeanFreePath(trk,prevstep,cond);
}

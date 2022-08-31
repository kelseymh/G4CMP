/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPQPDownconversion.hh
/// \brief Definition of the G4CMPQPDownconversion base class
// Converts Incoming Phonon to QPs

#include "G4CMPQPDownconversion.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPBogoliubov.hh"
#include "G4CMPSecondaryUtils.hh"
#include "G4CMPTrackUtils.hh"
#include "G4CMPUtils.hh"
#include "G4ExceptionSeverity.hh"
#include "G4LatticePhysical.hh"
#include "G4PhysicalConstants.hh"
#include "G4RandomDirection.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4VParticleChange.hh"
#include "G4CMPKaplanUtils.hh"
#include "Randomize.hh"
#include <cmath>

G4CMPQPDownconversion::G4CMPQPDownconversion(const G4String& aName)
  : G4VPhononProcess(aName, fPhononDownconversion)
{
}



G4CMPQPDownconversion::~G4CMPQPDownconversion() {;}

G4VParticleChange*
G4CMPQPDownconversion::PostStepDoIt(const G4Track& aTrack,
				    const G4Step& aStep) {
  //G4CMPBoundaryUtils::SetVerboseLevel(verboseLevel);

  aParticleChange.Initialize(aTrack);

  //Make sure phonon
  //if(!IsApplicable(aTrack.GetParticleDefinition)) {return &aParticleChange;}

  //Check if phonon hits superconducting film
  G4StepPoint* postStepPoint = aStep.GetPostStepPoint();
  if (postStepPoint->GetStepStatus() == fGeomBoundary ||
      postStepPoint->GetStepStatus() == fWorldBoundary) {
    return &aParticleChange;			// Don't want to reset IL
  }
  //Check if energy is sufficient to break QP
  const G4double phonon_E = aTrack.GetKineticEnergy();

  if(IsSubgap(phonon_E)) {
    //Do subgap things
  }

  //If yes, then (mayhaps) break
  MakeQPSecondaries(aTrack);
  

  //if downconverted, kill track
  if (aParticleChange.GetNumberOfSecondaries() > 0) {
    aParticleChange.ProposeEnergy(0.);
    aParticleChange.ProposeTrackStatus(fStopAndKill);

  //sanity check later
  }

  return &aParticleChange;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4CMPQPDownconversion::IsApplicable(const G4ParticleDefinition& aPD) {
  return 1;
}

void G4CMPQPDownconversion::MakeQPSecondaries(const G4Track& aTrack) {
 //Calc Energy + mom
 

  //Create 2 Tracks
  //G4Track* sec1 = G4CMP::CreatePhonon(aTrack.GetTouchable(), mode1,
				      //dir1, Esec1, aTrack.GetGlobalTime(),
                                      //aTrack.GetPosition());
  //G4Track* sec2 = G4CMP::CreateQP
  
  //Finalize change

  //aParticleChange.SetNumberOfSecondaries(2);
  //aParticleChange.AddSecondary(sec2);
  //aParticleChange.AddSecondary(sec1); 
}

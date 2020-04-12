/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// 20200331  G4CMP-196: Added impact ionization process

#include "G4CMPDriftImpactProcess.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPSecondaryUtils.hh"
#include "G4CMPUtils.hh"
#include "G4RandomDirection.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include <vector>


// Constructor and destructor

G4CMPDriftImpactProcess::
G4CMPDriftImpactProcess(const G4String &name)
  : G4CMPVDriftProcess(name, fImpactIonization) {}

G4CMPDriftImpactProcess::~G4CMPDriftImpactProcess() {
}


// Process actions

G4double 
G4CMPDriftImpactProcess::GetMeanFreePath(const G4Track&, G4double,
						G4ForceCondition* cond) {
  return  (IsElectron() ? G4CMPConfigManager::GetImpactLengthElectrons()
	   : IsHole()   ? G4CMPConfigManager::GetImpactLengthHoles()
	   : DBL_MAX);
}

G4VParticleChange* 
G4CMPDriftImpactProcess::PostStepDoIt(const G4Track& aTrack,
					     const G4Step& aStep) {
  aParticleChange.Initialize(aTrack);

  if (verboseLevel > 1) {
    G4cout << "G4CMPDriftImpactProcess::PostStepDoIt: "
           << aTrack.GetDefinition()->GetParticleName()
           << " impact ionization of an impurity."
           << G4endl;
  }

  // NOTE: If secondary gets real energy/momentum transfer, we *MUST*
  //       do the proper kinematics and update aParticleChange.

  // Create secondary with minimal energy, assuming no momentum transfer
  G4Track* knockon = G4CMP::CreateSecondary(aTrack, aTrack.GetDefinition(),
					    G4RandomDirection(), 1e-3*eV);
  aParticleChange.AddSecondary(knockon);

  ClearNumberOfInteractionLengthLeft();		// All processes should do this!
  return &aParticleChange;
}

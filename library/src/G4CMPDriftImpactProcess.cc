/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// 20200331  G4CMP-195/196: Added impact ionization and trapping

#include "G4CMPDriftImpactProcess.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4CMPSecondaryUtils.hh"
#include "G4CMPUtils.hh"
#include "G4LatticePhysical.hh"
#include "G4RandomDirection.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4Track.hh"
#include <vector>


// Constructor and destructor

G4CMPDriftImpactProcess::
G4CMPDriftImpactProcess(const G4String &name)
  : G4CMPVDriftProcess(name, fImpact) {}

G4CMPDriftImpactProcess::~G4CMPDriftImpactProcess() {
}


// Process actions

G4double 
G4CMPDriftImpactProcess::GetMeanFreePath(const G4Track&, G4double,
						G4ForceCondition* cond) {
  // CS TODO: How to access impactLength from ConfigManager?
  *cond = Forced;
  return DBL_MAX;
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

  // Create a like particle with random direction and 0 energy
  G4PrimaryParticle* primary = new G4PrimaryParticle();
  if (G4CMP::IsElectron(aTrack)) {
    primary->SetParticleDefinition(G4CMPDriftElectron::Definition());
  } else if (G4CMP::IsHole(aTrack)) {
    primary->SetParticleDefinition(G4CMPDriftHole::Definition());    
  } else {
    G4Exception("G4CMPLukeScattering::PostStepDoIt", "Luke002",
                EventMustBeAborted, "Unknown charge carrier");
    return &aParticleChange;
  }
  primary->SetMomentumDirection(G4RandomDirection());
  primary->SetKineticEnergy(0.);

  // Create new vertex
  G4PrimaryVertex* vertex;// = new G4PrimaryVertex(pos, time); // CS TODO: How do I access the pos and time of impact?
  vertex->SetPrimary(primary);

  // Add vertext to event
  //event->AddPrimaryVertex(vertex); // CS TODO: event is undefined. How do I add this to the current event?

  ClearNumberOfInteractionLengthLeft();		// All processes should do this!
  return &aParticleChange;
}

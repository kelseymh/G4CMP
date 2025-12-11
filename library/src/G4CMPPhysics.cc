/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
//
// Create particles and physics processes for phonons and charge carriers
// Usage:  [physics-list]->AddPhysics(new G4CMPPhysics);
//
// 20160901  M. Kelsey -- Add minimum-energy cut process
// 20161219  M. Kelsey -- Use particle table iterator directly
// 20170817  M. Kelsey -- Get verbosity from configuration
// 20170822  M. Kelsey -- Rename EnergyLimiter to TrackLimiter
// 20191017  M. Kelsey -- Add GenericIon to support energy partitioner
// 20200331  C. Stanford (G4CMP-195): Add charge trapping process
// 20200331  G4CMP-196: Added impact ionization process
// 20200426  G4CMP-196: Change "impact" to "trap ionization", separate
//		process instances for each beam/trap type.
// 20210203  G4CMP-241: SecondaryProduction must be last PostStep process.
// 20220331  G4CMP-293: Replace RegisterProcess() with local AddG4CMPProcess().
// 20251203  G4CMP-552: For G4 11, add hypernuclei daughters.
// 20251211  G4CMP-553: Add missing hypernuclei daughters (neutrinos).

#include "G4CMPPhysics.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftBoundaryProcess.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4CMPDriftRecombinationProcess.hh"
#include "G4CMPDriftTrappingProcess.hh"
#include "G4CMPDriftTrapIonization.hh"
#include "G4CMPInterValleyScattering.hh"
#include "G4CMPLukeScattering.hh"
#include "G4CMPPhononBoundaryProcess.hh"
#include "G4CMPSecondaryProduction.hh"
#include "G4CMPTimeStepper.hh"
#include "G4CMPTrackLimiter.hh"
#include "G4AntiLambda.hh"
#include "G4AntiNeutrinoE.hh"
#include "G4AntiNeutrinoMu.hh"
#include "G4AntiNeutron.hh"
#include "G4AntiProton.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4GenericIon.hh"
#include "G4Lambda.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4NeutrinoE.hh"
#include "G4NeutrinoMu.hh"
#include "G4Neutron.hh"
#include "G4ParticleTable.hh"
#include "G4PhononDownconversion.hh"
#include "G4PhononLong.hh"
#include "G4PhononScattering.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononTransSlow.hh"
#include "G4Positron.hh"
#include "G4PionMinus.hh"
#include "G4PionPlus.hh"
#include "G4PionZero.hh"
#include "G4ProcessManager.hh"
#include "G4Proton.hh"


// Constructor sets global verbosity

G4CMPPhysics::G4CMPPhysics(const G4String& name)
  : G4VPhysicsConstructor(name) {
  SetVerboseLevel(G4CMPConfigManager::GetVerboseLevel());
}

// Create phonon and charage carrier particles for later use

void G4CMPPhysics::ConstructParticle() {
  G4CMPDriftElectron::Definition();
  G4CMPDriftHole::Definition();
  G4PhononLong::Definition();
  G4PhononTransFast::Definition();
  G4PhononTransSlow::Definition();
  G4GenericIon::Definition();

  // For Geant4 11, define daughters of hypernuclei, created by GenericIon
  G4AntiLambda::Definition();
  G4AntiNeutrinoE::Definition();
  G4AntiNeutrinoMu::Definition();
  G4AntiNeutron::Definition();
  G4AntiProton::Definition();
  G4Electron::Definition();
  G4Gamma::Definition();
  G4Lambda::Definition();
  G4MuonMinus::Definition();
  G4MuonPlus::Definition();
  G4NeutrinoE::Definition();
  G4NeutrinoMu::Definition();
  G4Neutron::Definition();
  G4PionMinus::Definition();
  G4PionPlus::Definition();
  G4PionZero::Definition();
  G4Positron::Definition();
  G4Proton::Definition();
}

// Add physics processes to appropriate particles

void G4CMPPhysics::ConstructProcess() {
  // Only make processes once; will be deleted when physics list goes away
  G4VProcess* phScat  = new G4PhononScattering;
  G4VProcess* phRefl  = new G4CMPPhononBoundaryProcess;
  G4VProcess* phDown  = new G4PhononDownconversion;
  G4VProcess* tmStep  = new G4CMPTimeStepper;
  G4VProcess* driftB  = new G4CMPDriftBoundaryProcess;
  G4VProcess* ivScat  = new G4CMPInterValleyScattering;
  G4VProcess* luke    = new G4CMPLukeScattering(tmStep);
  G4VProcess* recomb  = new G4CMPDriftRecombinationProcess;
  G4VProcess* eLimit  = new G4CMPTrackLimiter;
  G4VProcess* trapping = new G4CMPDriftTrappingProcess;

  // NOTE: Trap ionization needs separate instances for each particle type
  G4ParticleDefinition* edrift = G4CMPDriftElectron::Definition();
  G4ParticleDefinition* hdrift = G4CMPDriftHole::Definition();

  G4VProcess* eeTrpI = new G4CMPDriftTrapIonization(edrift, edrift);
  G4VProcess* ehTrpI = new G4CMPDriftTrapIonization(edrift, hdrift);
  G4VProcess* heTrpI = new G4CMPDriftTrapIonization(hdrift, edrift);
  G4VProcess* hhTrpI = new G4CMPDriftTrapIonization(hdrift, hdrift);

  // Add processes only to locally known particles
  G4ParticleDefinition* particle = 0;

  particle = G4PhononLong::PhononDefinition();
  AddG4CMPProcess(phScat, particle);
  AddG4CMPProcess(phDown, particle);
  AddG4CMPProcess(phRefl, particle);
  AddG4CMPProcess(eLimit, particle);

  particle = G4PhononTransSlow::PhononDefinition();
  AddG4CMPProcess(phScat, particle);
  AddG4CMPProcess(phDown, particle);
  AddG4CMPProcess(phRefl, particle);
  AddG4CMPProcess(eLimit, particle);

  particle = G4PhononTransFast::PhononDefinition();
  AddG4CMPProcess(phScat, particle);
  AddG4CMPProcess(phDown, particle);
  AddG4CMPProcess(phRefl, particle);
  AddG4CMPProcess(eLimit, particle);

  particle = edrift;
  AddG4CMPProcess(tmStep, particle);
  AddG4CMPProcess(luke, particle);
  AddG4CMPProcess(ivScat, particle);
  AddG4CMPProcess(driftB, particle);
  AddG4CMPProcess(recomb, particle);
  AddG4CMPProcess(eLimit, particle);
  AddG4CMPProcess(trapping, particle);
  AddG4CMPProcess(eeTrpI, particle);	// e- projectile on both traps
  AddG4CMPProcess(ehTrpI, particle);

  particle = hdrift;
  AddG4CMPProcess(tmStep, particle);
  AddG4CMPProcess(luke, particle);
  AddG4CMPProcess(driftB, particle);
  AddG4CMPProcess(recomb, particle);
  AddG4CMPProcess(eLimit, particle);
  AddG4CMPProcess(trapping, particle);
  AddG4CMPProcess(heTrpI, particle);	// h+ projectile on both traps
  AddG4CMPProcess(hhTrpI, particle);

  AddSecondaryProduction();
}


// Register G4CMP processes here, instead of using G4CMPOrdParamTable.txt

void G4CMPPhysics::AddG4CMPProcess(G4VProcess* proc, G4ParticleDefinition* pd) {
  // Set process verbosity to match physics list, for diagnostics
  proc->SetVerboseLevel(verboseLevel);
  pd->GetProcessManager()->AddDiscreteProcess(proc, ordLast);
}


// Add charge and phonon generator to all charged particles

void G4CMPPhysics::AddSecondaryProduction() {
  G4VProcess* maker = new G4CMPSecondaryProduction;

  auto pIter = G4ParticleTable::GetParticleTable()->GetIterator();
  pIter->reset();
  while ((*pIter)()) {
    G4ParticleDefinition* particle = pIter->value();
    if (maker->IsApplicable(*particle)) AddG4CMPProcess(maker, particle); 
  }
}

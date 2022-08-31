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

#include "G4CMPPhysics.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftBoundaryProcess.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4CMPBogoliubov.hh"
#include "G4CMPDriftRecombinationProcess.hh"
#include "G4CMPDriftTrappingProcess.hh"
#include "G4CMPDriftTrapIonization.hh"
#include "G4CMPInterValleyScattering.hh"
#include "G4CMPLukeScattering.hh"
#include "G4CMPPhononBoundaryProcess.hh"
#include "G4CMPSecondaryProduction.hh"
#include "G4CMPTimeStepper.hh"
#include "G4CMPTrackLimiter.hh"
#include "G4GenericIon.hh"
#include "G4ParticleTable.hh"
#include "G4PhononDownconversion.hh"
#include "G4PhononLong.hh"
#include "G4PhononScattering.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononTransSlow.hh"
#include "G4CMPQPRecombination.hh"
#include "G4ProcessManager.hh"


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
  G4CMPBogoliubov::Definition();
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
  G4VProcess* qpRcmb = new G4CMPQPRecombination;

  // NOTE: Trap ionization needs separate instances for each particle type
  G4ParticleDefinition* edrift = G4CMPDriftElectron::Definition();
  G4ParticleDefinition* hdrift = G4CMPDriftHole::Definition();

  G4VProcess* eeTrpI = new G4CMPDriftTrapIonization(edrift, edrift);
  G4VProcess* ehTrpI = new G4CMPDriftTrapIonization(edrift, hdrift);
  G4VProcess* heTrpI = new G4CMPDriftTrapIonization(hdrift, edrift);
  G4VProcess* hhTrpI = new G4CMPDriftTrapIonization(hdrift, hdrift);

  // Set process verbosity to match physics list, for diagnostics
  if (verboseLevel>0) {
    phScat->SetVerboseLevel(verboseLevel);
    phRefl->SetVerboseLevel(verboseLevel);
    phDown->SetVerboseLevel(verboseLevel);
    tmStep->SetVerboseLevel(verboseLevel);
    driftB->SetVerboseLevel(verboseLevel);
    ivScat->SetVerboseLevel(verboseLevel);
    luke->SetVerboseLevel(verboseLevel);
    recomb->SetVerboseLevel(verboseLevel);
    eLimit->SetVerboseLevel(verboseLevel);
    trapping->SetVerboseLevel(verboseLevel);
    eeTrpI->SetVerboseLevel(verboseLevel);
    ehTrpI->SetVerboseLevel(verboseLevel);
    heTrpI->SetVerboseLevel(verboseLevel);
    hhTrpI->SetVerboseLevel(verboseLevel);
    qpRcmb->SetVerboseLevel(verboseLevel);
  }

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

  particle = G4CMPBogoliubov::Definition();
  RegisterProcess(qpRcmb, particle);

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

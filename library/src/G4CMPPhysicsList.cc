//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file exoticphysics/phonon/src/G4CMPPhysicsList.cc
/// \brief Implementation of the G4CMPPhysicsList class
//
// $Id$
//

#include "G4CMPPhysicsList.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4ios.hh"              

#include "G4PhononTransSlow.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononLong.hh"

#include "G4PhononScattering.hh"
#include "G4PhononDownconversion.hh"
#include "G4PhononReflection.hh"

#include "G4CMPDriftBoundaryProcess.hh"
#include "G4CMPhLukeScattering.hh"
#include "G4CMPeLukeScattering.hh"
// #include "EFieldProcess.hh"
// #include "ObliqueEFieldProcess.hh"
#include "G4CMPTimeStepper.hh"
#include "G4CMPInterValleyScattering.hh"

#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4StepLimiter.hh"
#include "G4UserSpecialCuts.hh"

#include "G4UserLimits.hh"



G4CMPPhysicsList::G4CMPPhysicsList(G4int verbose) :  G4VUserPhysicsList()
{
  if (verbose) G4cout << "G4CMPPhysicsList::constructor" << G4endl;

  defaultCutValue = DBL_MIN;	// 100*mm; 
  SetVerboseLevel(verbose);
}

G4CMPPhysicsList::~G4CMPPhysicsList() {}

void G4CMPPhysicsList::ConstructParticle() {
  G4CMPDriftElectron::Definition();
  G4CMPDriftHole::Definition();
  G4PhononLong::PhononDefinition();
  G4PhononTransFast::PhononDefinition();
  G4PhononTransSlow::PhononDefinition();
}


void G4CMPPhysicsList::ConstructProcess() {
  AddTransportation();

  // Only make common processes once (charge-carriers need separate instances)
  G4VProcess* phScat = new G4PhononScattering;
  G4VProcess* phRefl = new G4PhononReflection;
  G4VProcess* phDown = new G4PhononDownconversion;
  G4VProcess* tmStep = new G4CMPTimeStepper;
  G4VProcess* driftB = new G4CMPDriftBoundaryProcess;
  G4VProcess* ivScat = new G4CMPInterValleyScattering;

  // Set process verbosity to match physics list, for diagnostics
  phScat->SetVerboseLevel(verboseLevel);
  phRefl->SetVerboseLevel(verboseLevel);
  phDown->SetVerboseLevel(verboseLevel);
  tmStep->SetVerboseLevel(verboseLevel);
  driftB->SetVerboseLevel(verboseLevel);
  ivScat->SetVerboseLevel(verboseLevel);

  G4ParticleDefinition* particle = 0;
  G4ProcessManager* pmanager = 0;

  // Add processes only to locally known particles
  particle = G4PhononLong::PhononDefinition();
  pmanager = particle->GetProcessManager();
  pmanager->AddDiscreteProcess(phScat);
  pmanager->AddDiscreteProcess(phDown);
  pmanager->AddDiscreteProcess(phRefl);

  particle = G4PhononTransSlow::PhononDefinition();
  pmanager = particle->GetProcessManager();
  pmanager->AddDiscreteProcess(phScat);
  pmanager->AddDiscreteProcess(phDown);
  pmanager->AddDiscreteProcess(phRefl);

  particle = G4PhononTransFast::PhononDefinition();
  pmanager = particle->GetProcessManager();
  pmanager->AddDiscreteProcess(phScat);
  pmanager->AddDiscreteProcess(phDown);
  pmanager->AddDiscreteProcess(phRefl);

  particle = G4CMPDriftElectron::Definition();
  pmanager = particle->GetProcessManager();
  pmanager->AddDiscreteProcess(tmStep);
  pmanager->AddDiscreteProcess(new G4CMPeLukeScattering(tmStep));
  pmanager->AddDiscreteProcess(ivScat);
  pmanager->AddDiscreteProcess(driftB);

  particle = G4CMPDriftHole::Definition();
  pmanager = particle->GetProcessManager();
  pmanager->AddDiscreteProcess(tmStep);
  pmanager->AddDiscreteProcess(new G4CMPhLukeScattering(tmStep));
  pmanager->AddDiscreteProcess(driftB);
}

void G4CMPPhysicsList::SetCuts()
{
  // These values are used as the default production thresholds
  // for the world volume.
  SetCutsWithDefault();
}



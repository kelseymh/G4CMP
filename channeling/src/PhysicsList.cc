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
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// General
#include "PhysicsList.hh"
#include "globals.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ProcessManager.hh"

// Single Scattering EM model
#include "PhysListEmStandardSS.hh"

// Decay
#include "G4Decay.hh"
#include "StepMax.hh"


// Particle
#include "G4ParticleDefinition.hh"

// Bosons
#include "G4ChargedGeantino.hh"
#include "G4Geantino.hh"
#include "G4Gamma.hh"
#include "G4OpticalPhoton.hh"

// Leptons
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4NeutrinoMu.hh"
#include "G4AntiNeutrinoMu.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoE.hh"

// Hadrons
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"

// Channeling
#include "ProcessChanneling.hh"

PhysicsList::PhysicsList():  G4VUserPhysicsList(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess(){
    AddTransportation();
    //AddDecay();
    AddStepMax();
    AddChanneling();
    //AddStandardSS();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle(){
    // pseudo-particles
    G4Geantino::GeantinoDefinition();
    G4ChargedGeantino::ChargedGeantinoDefinition();
    
    // gamma
    G4Gamma::GammaDefinition();
    
    // leptons
    G4Electron::ElectronDefinition();
    G4Positron::PositronDefinition();
    G4MuonPlus::MuonPlusDefinition();
    G4MuonMinus::MuonMinusDefinition();
    
    G4NeutrinoE::NeutrinoEDefinition();
    G4AntiNeutrinoE::AntiNeutrinoEDefinition();
    G4NeutrinoMu::NeutrinoMuDefinition();
    G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
    
    // mesons
    G4MesonConstructor mConstructor;
    mConstructor.ConstructParticle();
    
    // barions
    G4BaryonConstructor bConstructor;
    bConstructor.ConstructParticle();
    
    // ions
    G4IonConstructor iConstructor;
    iConstructor.ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddStandardSS(){
    G4VPhysicsConstructor*  fEmPhysicsList = new PhysListEmStandardSS("standardSS");
    fEmPhysicsList->ConstructProcess();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddChanneling(){
    theParticleIterator->reset();
    while( (*theParticleIterator)() ){
        G4ParticleDefinition* particle = theParticleIterator->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();
        G4String particleName = particle->GetParticleName();
        
        if(particleName=="proton")
        {
            G4cout<<"\n\nPhysicsList::ConstructParticle: \n\tFOUND PROTON...\n"<<G4endl;
            pmanager->AddDiscreteProcess(new ProcessChanneling());
        }
    }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddDecay(){
    // Add Decay Process
    G4Decay* fDecayProcess = new G4Decay();
    
    theParticleIterator->reset();
    while( (*theParticleIterator)() ){
        G4ParticleDefinition* particle = theParticleIterator->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();
        
        if (fDecayProcess->IsApplicable(*particle) && !particle->IsShortLived()) {
            
            pmanager->AddProcess(fDecayProcess);
            
            // set ordering for PostStepDoIt and AtRestDoIt
            pmanager->SetProcessOrdering(fDecayProcess, idxPostStep);
            pmanager->SetProcessOrdering(fDecayProcess, idxAtRest);
            
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddStepMax()
{
    // Step limitation seen as a process
    StepMax* stepMaxProcess = new StepMax();
    
    theParticleIterator->reset();
    while ((*theParticleIterator)()){
        G4ParticleDefinition* particle = theParticleIterator->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();
        
        if (stepMaxProcess->IsApplicable(*particle) && !particle->IsShortLived())
        {
            pmanager ->AddDiscreteProcess(stepMaxProcess);
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void PhysicsList::SetCuts()
{
    // These values are used as the default production thresholds
    // for the world volume.
    SetCutsWithDefault();
    
}



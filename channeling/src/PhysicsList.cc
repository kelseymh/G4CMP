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

#include "PhysicsList.hh"

// General
#include "globals.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessVector.hh"

// Particle
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

// Single Scattering EM model
#include "G4CoulombScattering.hh"

// Decay
#include "G4Decay.hh"
#include "StepMax.hh"

// Channeling
#include "ProcessChanneling.hh"

#include "XVCrystalCharacteristic.hh"
#include "XVCrystalIntegratedDensity.hh"

#include "XCrystalPlanarMoliereTempPotential.hh"
#include "XCrystalPlanarMoliereElectricField.hh"
#include "XCrystalPlanarNucleiDensity.hh"
#include "XCrystalPlanarMoliereElectronDensity.hh"
#include "XCrystalCharacteristicArray.hh"
#include "XCrystalIntegratedDensityPlanar.hh"

// Proton hadronic
#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4ProtonInelasticProcess.hh"

#include "G4LElastic.hh"
#include "G4LFission.hh"
#include "G4LCapture.hh"

#include "G4TheoFSGenerator.hh"
#include "G4ExcitationHandler.hh"
#include "G4Evaporation.hh"
#include "G4CompetitiveFission.hh"
#include "G4FermiBreakUp.hh"
#include "G4StatMF.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4Fancy3DNucleus.hh"
#include "G4LEProtonInelastic.hh"
#include "G4StringModel.hh"
#include "G4PreCompoundModel.hh"
#include "G4FTFModel.hh"
#include "G4QGSMFragmentation.hh"
#include "G4ExcitedStringDecay.hh"

#include "G4CHIPSElastic.hh"
#include "G4CascadeInterface.hh"

#include "G4TheoFSGenerator.hh"
#include "G4ExcitationHandler.hh"
#include "G4PreCompoundModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4FTFModel.hh"
#include "G4LundStringFragmentation.hh"
#include "G4ExcitedStringDecay.hh"

#include "G4CHIPSElasticXS.hh"
#include "G4BGGNucleonInelasticXS.hh"
#include "G4ProtonInelasticCrossSection.hh"
#include "G4CrossSectionPairGG.hh"

// Wrapper
#include "XWrapperProcess.hh"

PhysicsList::PhysicsList():  G4VUserPhysicsList(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess(){
    AddTransportation();
    //AddChanneling();
    AddInelaticProcesses();
    AddStandardSS();
    AddDecay();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle(){
    G4BosonConstructor pBosonConstructor;
    pBosonConstructor.ConstructParticle();
    
    G4LeptonConstructor pLeptonConstructor;
    pLeptonConstructor.ConstructParticle();
    
    G4MesonConstructor pMesonConstructor;
    pMesonConstructor.ConstructParticle();
    
    G4BaryonConstructor pBaryonConstructor;
    pBaryonConstructor.ConstructParticle();
    
    G4IonConstructor pIonConstructor;
    pIonConstructor.ConstructParticle();
    
    G4ShortLivedConstructor pShortLivedConstructor;
    pShortLivedConstructor.ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PhysicsList::AddInelaticProcesses(){
    // this will be the model class for high energies
    G4TheoFSGenerator * theTheoModel = new G4TheoFSGenerator;
    
    // all models for treatment of thermal nucleus
    G4Evaporation * theEvaporation = new G4Evaporation;
    G4FermiBreakUp * theFermiBreakUp = new G4FermiBreakUp;
    G4StatMF * theMF = new G4StatMF;
    
    // Evaporation logic
    G4ExcitationHandler * theHandler = new G4ExcitationHandler;
    theHandler->SetEvaporation(theEvaporation);
    theHandler->SetFermiModel(theFermiBreakUp);
    theHandler->SetMultiFragmentation(theMF);
    theHandler->SetMaxAandZForFermiBreakUp(12, 6);
    theHandler->SetMinEForMultiFrag(3*MeV);
    
    // Pre equilibrium stage
    G4PreCompoundModel * thePreEquilib = new G4PreCompoundModel(theHandler);
    
    // a no-cascade generator-precompound interaface
    G4GeneratorPrecompoundInterface * theCascade = new G4GeneratorPrecompoundInterface;
    theCascade->SetDeExcitation(thePreEquilib);
    
    // here come the high energy parts
    // the string model; still not quite according to design - Explicite use of the forseen interfaces
    // will be tested and documented in this program by beta-02 at latest.
    G4VPartonStringModel * theStringModel;
    theStringModel = new G4FTFModel;
    theTheoModel->SetTransport(theCascade);
    theTheoModel->SetHighEnergyGenerator(theStringModel);
    theTheoModel->SetMinEnergy(19*GeV);
    theTheoModel->SetMaxEnergy(100*TeV);
    
    G4VLongitudinalStringDecay * theFragmentation = new G4QGSMFragmentation;
    G4ExcitedStringDecay * theStringDecay = new G4ExcitedStringDecay(theFragmentation);
    theStringModel->SetFragmentationModel(theStringDecay);
    
    // done with the generator model (most of the above is also available as default)
    G4HadronElasticProcess* theElasticProcess =
    new G4HadronElasticProcess;
    G4LElastic* theElasticModel = new G4LElastic;
    theElasticProcess->RegisterMe(theElasticModel);
    
    theParticleIterator->reset();
    while( (*theParticleIterator)() ){
        G4ParticleDefinition* particle = theParticleIterator->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();
        G4String particleName = particle->GetParticleName();
        if (particleName == "proton" ) {
            G4ProtonInelasticProcess* theInelasticProcess =
            new G4ProtonInelasticProcess("inelastic");
            G4LEProtonInelastic* theInelasticModel = new G4LEProtonInelastic;
            theInelasticProcess->RegisterMe(theInelasticModel);
            theInelasticProcess->RegisterMe(theTheoModel);

            XWrapperProcess *elasticProcess_wrapper = new XWrapperProcess();
            elasticProcess_wrapper->registerProcess(theElasticProcess);
            XWrapperProcess *inelasticProcess_wrapper = new XWrapperProcess();
            inelasticProcess_wrapper->registerProcess(theInelasticProcess);
            
            //pmanager->AddDiscreteProcess(theElasticProcess);
            //pmanager->AddDiscreteProcess(theInelasticProcess);
            pmanager->AddDiscreteProcess(elasticProcess_wrapper);
            pmanager->AddDiscreteProcess(inelasticProcess_wrapper);
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddStandardSS(){
    
    theParticleIterator->reset();
    while( (*theParticleIterator)() ){
        G4ParticleDefinition* particle = theParticleIterator->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();
        G4String particleName = particle->GetParticleName();
        if (particleName == "proton" ) {
            G4CoulombScattering* cs = new G4CoulombScattering();
            cs->SetBuildTableFlag(false);
            
            XWrapperProcess *cs_wrapper = new XWrapperProcess();
            cs_wrapper->registerProcess(cs);

            //pmanager->AddDiscreteProcess(cs);
            pmanager->AddDiscreteProcess(cs_wrapper);
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddChanneling(){

    XVCrystalCharacteristic* vPotentialEnergy = new XCrystalPlanarMoliereTempPotential();
    
    XVCrystalCharacteristic* vElectricField = new XCrystalPlanarMoliereElectricField();
    
    XVCrystalCharacteristic* vNucleiDensity = new XCrystalPlanarNucleiDensity();
    XVCrystalCharacteristic* vElectronDensity = new XCrystalPlanarMoliereElectronDensity();
    XVCrystalIntegratedDensity* vIntegratedDensity = new XCrystalIntegratedDensityPlanar();
    vIntegratedDensity->SetPotential(vPotentialEnergy);
    vIntegratedDensity->SetDensity(vNucleiDensity);
    
    theParticleIterator->reset();
    while( (*theParticleIterator)() ){
        G4ParticleDefinition* particle = theParticleIterator->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();
        G4String particleName = particle->GetParticleName();
        
        if(particleName == "proton")
        {
            G4cout<<"\n\nPhysicsList::ConstructParticle: \n\tFOUND PROTON...\n"<<G4endl;
            ProcessChanneling* channeling =  new ProcessChanneling();
            channeling->SetPotential(vPotentialEnergy);
            channeling->SetIntegratedDensity(vIntegratedDensity);
            
            pmanager->AddDiscreteProcess(channeling);
            
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



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
#include "G4RunManager.hh"

// Particle
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

// EM
#include "G4CoulombScattering.hh"
#include "G4hIonisation.hh"
#include "G4ScreenedNuclearRecoil.hh"

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
#include "XWrapperDiscreteProcess.hh"
#include "XWrapperContinuousDiscreteProcess.hh"

PhysicsList::PhysicsList():  G4VUserPhysicsList(){
    fFileName = "";
    
    fScatteringType = "ss";
    
    bChannelingOn = true;
    
    bWrapperOn = true;
    
    fMessenger = new PhysicsListMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess(){
    AddTransportation();
    
    if(bChannelingOn){
        AddChanneling();
    }
    
    AddInelaticProcesses();
    
    AddDecay();

    if(fScatteringType == "ss"){
        AddStandardSS();
    }
    else if(fScatteringType == "nr"){
        AddStandardNR();
    }
    else{
        G4cout << "PhysicsList::ConstructProcess - No Scattering model selected" << G4endl;
    }
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
    G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
    G4LElastic* theElasticModel = new G4LElastic;
    theElasticProcess->RegisterMe(theElasticModel);
    
    theParticleIterator->reset();
    while( (*theParticleIterator)() ){
        G4ParticleDefinition* particle = theParticleIterator->value();
        G4ProcessManager* pManager = particle->GetProcessManager();
        G4String particleName = particle->GetParticleName();
        if (particleName == "proton" ) {
            G4ProtonInelasticProcess* theInelasticProcess = new G4ProtonInelasticProcess("inelastic");
            G4LEProtonInelastic* theInelasticModel = new G4LEProtonInelastic;
            theInelasticProcess->RegisterMe(theInelasticModel);
            theInelasticProcess->RegisterMe(theTheoModel);
            
            XWrapperDiscreteProcess *elasticProcess_wrapper = new XWrapperDiscreteProcess();
            elasticProcess_wrapper->RegisterProcess(theElasticProcess);
            XWrapperDiscreteProcess *inelasticProcess_wrapper = new XWrapperDiscreteProcess();
            inelasticProcess_wrapper->RegisterProcess(theInelasticProcess);
            
            if(bWrapperOn){
                pManager->AddDiscreteProcess(elasticProcess_wrapper,1);
                pManager->AddDiscreteProcess(inelasticProcess_wrapper,1);
                G4cout<<"\nPhysicsList::AddInelaticProcesses: Hadronic inelastic process WITH WRAPPER added...\n"<<G4endl;
                
            }
            else{
                pManager->AddDiscreteProcess(theElasticProcess);
                pManager->AddDiscreteProcess(theInelasticProcess);
                G4cout<<"\nPhysicsList::AddInelaticProcesses: Hadronic inelastic process WITHOUT WRAPPER added...\n"<<G4endl;
                
            }
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddStandardSS(){
    
    theParticleIterator->reset();
    while( (*theParticleIterator)() ){
        G4ParticleDefinition* particle = theParticleIterator->value();
        G4ProcessManager* pManager = particle->GetProcessManager();
        G4String particleName = particle->GetParticleName();
        if (particleName == "proton" ) {
            
            G4CoulombScattering* cs = new G4CoulombScattering();
            cs->SetBuildTableFlag(false);
            XWrapperDiscreteProcess *cs_wrapper = new XWrapperDiscreteProcess();
            cs_wrapper->RegisterProcess(cs,1);
            
            G4hIonisation* theppIonisation = new G4hIonisation();
            theppIonisation->SetStepFunction(0.05, 1*um);
            XWrapperContinuousDiscreteProcess *hionisation_wrapper = new XWrapperContinuousDiscreteProcess();
            hionisation_wrapper->RegisterProcess(theppIonisation,-1);
            
            if(bWrapperOn){
                pManager->AddDiscreteProcess(cs_wrapper,1);
                pManager->AddProcess(hionisation_wrapper,-1, 2, 2);
                
            }
            else{
                
                pManager->AddDiscreteProcess(cs,1);
                pManager->AddProcess(theppIonisation,-1, 1, 1);
                G4cout<<"\nPhysicsList::AddStandardSS: Single Scattering process WITHOUT WRAPPER added...\n"<<G4endl;
            }
            
            G4cout<< "\nPhysicsList::AddStandardSS: Single Scattering process added";
            if(bWrapperOn){
                G4cout << "WITH WRAPPER...\n"<<G4endl;
            }
            else{
                G4cout << "WITHOUT WRAPPER...\n"<<G4endl;
            }
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddStandardNR(){
    
    theParticleIterator->reset();
    while( (*theParticleIterator)() ){
        G4ParticleDefinition* particle = theParticleIterator->value();
        G4ProcessManager* pManager = particle->GetProcessManager();
        G4String particleName = particle->GetParticleName();
        if (particleName == "proton" ) {
            
            G4ScreenedNuclearRecoil* nucr = new G4ScreenedNuclearRecoil();
            XWrapperDiscreteProcess *nucr_wrapper = new XWrapperDiscreteProcess();
            nucr_wrapper->RegisterProcess(nucr,1);
            
            G4hIonisation* theppIonisation = new G4hIonisation();
            theppIonisation->SetStepFunction(0.05, 1*um);
            XWrapperContinuousDiscreteProcess *hionisation_wrapper = new XWrapperContinuousDiscreteProcess();
            hionisation_wrapper->RegisterProcess(theppIonisation,-1);
            
            if(bWrapperOn){
                pManager->AddDiscreteProcess(nucr_wrapper,1);
                pManager->AddProcess(hionisation_wrapper,-1, 2, 2);
            }
            else{
                
                pManager->AddDiscreteProcess(nucr,1);
                pManager->AddProcess(theppIonisation,-1, 2, 2);
            }
            
            G4cout<< "\nPhysicsList::AddStandardNR: Nuclear Recoil Scattering process added";
            if(bWrapperOn){
                G4cout << "WITH WRAPPER...\n" << G4endl;
            }
            else{
                G4cout << "WITHOUT WRAPPER...\n" << G4endl;
            }
            
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddChanneling(){
    
    XVCrystalCharacteristic* vPotentialEnergy = new XCrystalPlanarMoliereTempPotential();
    
    XVCrystalCharacteristic* vElectricField = new XCrystalPlanarMoliereElectricField();
    
    XVCrystalCharacteristic* vNucleiDensity = new XCrystalPlanarNucleiDensity();
    XVCrystalCharacteristic* vElectronDensity = new XCrystalPlanarMoliereElectronDensity();
    
    XVCrystalIntegratedDensity* vIntegratedDensityNuclei = new XCrystalIntegratedDensityPlanar();
    vIntegratedDensityNuclei->SetPotential(vPotentialEnergy);
    vIntegratedDensityNuclei->SetDensity(vNucleiDensity);
    
    XVCrystalIntegratedDensity* vIntegratedDensityElectron = new XCrystalIntegratedDensityPlanar();
    vIntegratedDensityElectron->SetPotential(vPotentialEnergy);
    vIntegratedDensityElectron->SetDensity(vElectronDensity);
    
    theParticleIterator->reset();
    while( (*theParticleIterator)() ){
        G4ParticleDefinition* particle = theParticleIterator->value();
        G4ProcessManager* pManager = particle->GetProcessManager();
        G4String particleName = particle->GetParticleName();
        
        if(particleName == "proton")
        {
            ProcessChanneling* channeling =  new ProcessChanneling();
            channeling->SetPotential(vPotentialEnergy);
            channeling->SetIntegratedDensityNuclei(vIntegratedDensityNuclei);
            channeling->SetIntegratedDensityElectron(vIntegratedDensityElectron);
            channeling->SetFileName(fFileName);
            pManager->AddDiscreteProcess(channeling);
            G4cout<<"\nPhysicsList::AddChanneling: Channeling process added...\n"<<G4endl;
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
        G4ProcessManager* pManager = particle->GetProcessManager();
        
        if (fDecayProcess->IsApplicable(*particle) && !particle->IsShortLived()) {
            
            pManager->AddProcess(fDecayProcess);
            
            // set ordering for PostStepDoIt and AtRestDoIt
            pManager->SetProcessOrdering(fDecayProcess, idxPostStep);
            pManager->SetProcessOrdering(fDecayProcess, idxAtRest);
        }
    }
    G4cout<<"\nPhysicsList::AddDecay: Decay processes sdded...\n"<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddStepMax()
{
    // Step limitation seen as a process
    StepMax* stepMaxProcess = new StepMax();
    
    theParticleIterator->reset();
    while ((*theParticleIterator)()){
        G4ParticleDefinition* particle = theParticleIterator->value();
        G4ProcessManager* pManager = particle->GetProcessManager();
        
        if (stepMaxProcess->IsApplicable(*particle) && !particle->IsShortLived())
        {
            pManager ->AddDiscreteProcess(stepMaxProcess);
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::SetFileName(const G4String& vFilename){
    if(fFileName != vFilename){
        fFileName = vFilename;
        G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4String PhysicsList::GetFileName(){
    return fFileName;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::SetScatteringType(const G4String& vScatteringType){
    if(fScatteringType != vScatteringType){
        fScatteringType = vScatteringType;
        G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4String PhysicsList::GetScatteringType(){
    return fScatteringType;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::EnableChanneling(G4bool flag) {
    if(bChannelingOn != flag){
        bChannelingOn = flag;
        G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::EnableWrapper(G4bool flag) {
    if(bWrapperOn != flag){
        bWrapperOn = flag;
        G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


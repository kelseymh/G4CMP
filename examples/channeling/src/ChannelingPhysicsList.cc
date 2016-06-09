/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#include "ChannelingPhysicsList.hh"

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

#include "G4eCoulombScatteringModel.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4Generator2BS.hh"
#include "G4SeltzerBergerModel.hh"
#include "G4PenelopeIonisationModel.hh"
#include "G4UniversalFluctuation.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4RayleighScattering.hh"
#include "G4PEEffectFluoModel.hh"
#include "G4KleinNishinaModel.hh"
#include "G4LowEPComptonModel.hh"
#include "G4PenelopeGammaConversionModel.hh"
#include "G4LivermorePhotoElectricModel.hh"

#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4NuclearStopping.hh"

// Decay
#include "G4Decay.hh"
#include "ChannelingStepLimiter.hh"

// Channeling
#include "G4CMPChanneling.hh"
#include "G4CMPChannelingSlow.hh"

#include "XVCrystalCharacteristic.hh"
#include "XVCrystalIntegratedDensity.hh"

#include "XCrystalPlanarMoliereTempPotential.hh"
#include "XCrystalPlanarMoliereElectricField.hh"
#include "XCrystalPlanarNucleiDensity.hh"
#include "XCrystalPlanarMoliereElectronDensity.hh"
#include "XCrystalCharacteristicArray.hh"
#include "XCrystalIntegratedDensityPlanar.hh"
#include "XCrystalIntegratedDensityHub.hh"

// Hadronic
#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"

#include "G4ProtonInelasticProcess.hh"
#include "G4PionPlusInelasticProcess.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4KaonPlusInelasticProcess.hh"
#include "G4KaonZeroSInelasticProcess.hh"
#include "G4KaonZeroLInelasticProcess.hh"
#include "G4KaonMinusInelasticProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4AntiProtonInelasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4AntiNeutronInelasticProcess.hh"
#include "G4LambdaInelasticProcess.hh"
#include "G4AntiLambdaInelasticProcess.hh"
#include "G4SigmaPlusInelasticProcess.hh"
#include "G4SigmaMinusInelasticProcess.hh"
#include "G4AntiSigmaPlusInelasticProcess.hh"
#include "G4AntiSigmaMinusInelasticProcess.hh"
#include "G4XiZeroInelasticProcess.hh"
#include "G4XiMinusInelasticProcess.hh"
#include "G4AntiXiZeroInelasticProcess.hh"
#include "G4AntiXiMinusInelasticProcess.hh"
#include "G4DeuteronInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4AlphaInelasticProcess.hh"
#include "G4OmegaMinusInelasticProcess.hh"
#include "G4AntiOmegaMinusInelasticProcess.hh"

#include "G4HadronElastic.hh"
#include "G4LFission.hh"
#include "G4NeutronRadCapture.hh"

#include "G4LEPionPlusInelastic.hh"
#include "G4LEPionMinusInelastic.hh"
#include "G4LEKaonPlusInelastic.hh"
#include "G4LEKaonZeroSInelastic.hh"
#include "G4LEKaonZeroLInelastic.hh"
#include "G4LEKaonMinusInelastic.hh"
#include "G4LEProtonInelastic.hh"
#include "G4LEAntiProtonInelastic.hh"
#include "G4LENeutronInelastic.hh"
#include "G4LEAntiNeutronInelastic.hh"
#include "G4LELambdaInelastic.hh"
#include "G4LEAntiLambdaInelastic.hh"
#include "G4LESigmaPlusInelastic.hh"
#include "G4LESigmaMinusInelastic.hh"
#include "G4LEAntiSigmaPlusInelastic.hh"
#include "G4LEAntiSigmaMinusInelastic.hh"
#include "G4LEXiZeroInelastic.hh"
#include "G4LEXiMinusInelastic.hh"
#include "G4LEAntiXiZeroInelastic.hh"
#include "G4LEAntiXiMinusInelastic.hh"
#include "G4LEDeuteronInelastic.hh"
#include "G4LETritonInelastic.hh"
#include "G4LEAlphaInelastic.hh"
#include "G4LEOmegaMinusInelastic.hh"
#include "G4LEAntiOmegaMinusInelastic.hh"

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

#include "G4QMDReaction.hh"
#include "G4BinaryLightIonReaction.hh"
#include "G4CascadeInterface.hh"

#include "G4LundStringFragmentation.hh"

//#include "G4CHIPSElasticXS.hh"
#include "G4BGGNucleonInelasticXS.hh"
#include "G4ProtonInelasticCrossSection.hh"
#include "G4CrossSectionPairGG.hh"

// Wrapper
#include "ChannelingDiscreteWrapper.hh"
#include "ChannelingContinuousDiscreteWrapper.hh"

ChannelingPhysicsList::ChannelingPhysicsList():  G4VUserPhysicsList(){
    fFileName = "";
    
    fScatteringType = "ss";
    
    bChannelingOn = true;
    
    bWrapperOn = true;
    
    bDecayOn = true;
    
    fMessenger = new ChannelingPhysicsMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ChannelingPhysicsList::~ChannelingPhysicsList(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChannelingPhysicsList::ConstructProcess(){
    AddTransportation();
    
    if(bChannelingOn){
        AddChanneling();
    }
    
    if(bDecayOn){
        AddInelaticProcesses();
        AddDecay();
    }
        
    if(fScatteringType != "nosc"){
        AddEM();
    }
    
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChannelingPhysicsList::ConstructParticle(){
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
void ChannelingPhysicsList::AddInelaticProcesses(){
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
    G4HadronElastic* theElasticModel = new G4HadronElastic;
    theElasticProcess->RegisterMe(theElasticModel);
    
    ChannelingDiscreteWrapper *theElasticProcess_wrapper = new ChannelingDiscreteWrapper();
    theElasticProcess_wrapper->RegisterProcess(theElasticProcess);
    
    theParticleIterator->reset();
    while( (*theParticleIterator)() ){
        G4ParticleDefinition* particle = theParticleIterator->value();
        G4ProcessManager* pManager = particle->GetProcessManager();
        G4String particleName = particle->GetParticleName();
        if (particleName == "proton" ) {
            G4ProtonInelasticProcess* theInelasticProcess =
            new G4ProtonInelasticProcess("inelastic");
            G4LEProtonInelastic* theInelasticModel =
            new G4LEProtonInelastic;
            theInelasticProcess->RegisterMe(theInelasticModel);
            theInelasticProcess->RegisterMe(theTheoModel);
            
            ChannelingDiscreteWrapper *theInelasticProcess_wrapper = new ChannelingDiscreteWrapper();
            theInelasticProcess_wrapper->RegisterProcess(theInelasticProcess);
            
            if(bWrapperOn){
                pManager->AddDiscreteProcess(theElasticProcess_wrapper);
                pManager->AddDiscreteProcess(theInelasticProcess_wrapper);
            }
            else{
                pManager->AddDiscreteProcess(theElasticProcess);
                pManager->AddDiscreteProcess(theInelasticProcess);
            }
        }
        else if (particleName == "pi+") {
            G4PionPlusInelasticProcess* theInelasticProcess =
            new G4PionPlusInelasticProcess("inelastic");
            G4LEPionPlusInelastic* theInelasticModel =
            new G4LEPionPlusInelastic;
            theInelasticProcess->RegisterMe(theInelasticModel);
            theInelasticProcess->RegisterMe(theTheoModel);
            
            ChannelingDiscreteWrapper *theInelasticProcess_wrapper = new ChannelingDiscreteWrapper();
            theInelasticProcess_wrapper->RegisterProcess(theInelasticProcess);
            
            if(bWrapperOn){
                pManager->AddDiscreteProcess(theElasticProcess_wrapper);
                pManager->AddDiscreteProcess(theInelasticProcess_wrapper);
            }
            else{
                pManager->AddDiscreteProcess(theElasticProcess);
                pManager->AddDiscreteProcess(theInelasticProcess);
            }
        }
        else if (particleName == "pi-") {
            G4PionMinusInelasticProcess* theInelasticProcess =
            new G4PionMinusInelasticProcess("inelastic");
            G4LEPionMinusInelastic* theInelasticModel =
            new G4LEPionMinusInelastic;
            theInelasticProcess->RegisterMe(theInelasticModel);
            theInelasticProcess->RegisterMe(theTheoModel);
            
            ChannelingDiscreteWrapper *theInelasticProcess_wrapper = new ChannelingDiscreteWrapper();
            theInelasticProcess_wrapper->RegisterProcess(theInelasticProcess);
            
            if(bWrapperOn){
                pManager->AddDiscreteProcess(theElasticProcess_wrapper);
                pManager->AddDiscreteProcess(theInelasticProcess_wrapper);
            }
            else{
                pManager->AddDiscreteProcess(theElasticProcess);
                pManager->AddDiscreteProcess(theInelasticProcess);
            }
        }
        else if (particleName == "kaon+") {
            
            G4KaonPlusInelasticProcess* theInelasticProcess =
            new G4KaonPlusInelasticProcess("inelastic");
            G4LEKaonPlusInelastic* theInelasticModel = new G4LEKaonPlusInelastic;
            theInelasticProcess->RegisterMe(theInelasticModel);
            theInelasticProcess->RegisterMe(theTheoModel);
            
            ChannelingDiscreteWrapper *theInelasticProcess_wrapper = new ChannelingDiscreteWrapper();
            theInelasticProcess_wrapper->RegisterProcess(theInelasticProcess);
            
            if(bWrapperOn){
                pManager->AddDiscreteProcess(theElasticProcess_wrapper);
                pManager->AddDiscreteProcess(theInelasticProcess_wrapper);
            }
            else{
                pManager->AddDiscreteProcess(theElasticProcess);
                pManager->AddDiscreteProcess(theInelasticProcess);
            }
        }
        else if (particleName == "kaon0S") {
            
            G4KaonZeroSInelasticProcess* theInelasticProcess =
            new G4KaonZeroSInelasticProcess("inelastic");
            G4LEKaonZeroSInelastic* theInelasticModel =
            new G4LEKaonZeroSInelastic;
            theInelasticProcess->RegisterMe(theInelasticModel);
            theInelasticProcess->RegisterMe(theTheoModel);
            
            pManager->AddDiscreteProcess(theElasticProcess);
            pManager->AddDiscreteProcess(theInelasticProcess);
        }
        else if (particleName == "kaon0L") {
            
            G4KaonZeroLInelasticProcess* theInelasticProcess =
            new G4KaonZeroLInelasticProcess("inelastic");
            G4LEKaonZeroLInelastic* theInelasticModel =
            new G4LEKaonZeroLInelastic;
            theInelasticProcess->RegisterMe(theInelasticModel);
            theInelasticProcess->RegisterMe(theTheoModel);
            
            pManager->AddDiscreteProcess(theElasticProcess);
            pManager->AddDiscreteProcess(theInelasticProcess);
            
        }
        else if (particleName == "kaon-") {
            
            G4KaonMinusInelasticProcess* theInelasticProcess =
            new G4KaonMinusInelasticProcess("inelastic");
            G4LEKaonMinusInelastic* theInelasticModel =
            new G4LEKaonMinusInelastic;
            theInelasticProcess->RegisterMe(theInelasticModel);
            theInelasticProcess->RegisterMe(theTheoModel);
            
            ChannelingDiscreteWrapper *theInelasticProcess_wrapper = new ChannelingDiscreteWrapper();
            theInelasticProcess_wrapper->RegisterProcess(theInelasticProcess);
            
            if(bWrapperOn){
                pManager->AddDiscreteProcess(theElasticProcess_wrapper);
                pManager->AddDiscreteProcess(theInelasticProcess_wrapper);
            }
            else{
                pManager->AddDiscreteProcess(theElasticProcess);
                pManager->AddDiscreteProcess(theInelasticProcess);
            }
        }
        else if (particleName == "anti_proton") {
            
            G4AntiProtonInelasticProcess* theInelasticProcess =
            new G4AntiProtonInelasticProcess("inelastic");
            G4LEAntiProtonInelastic* theInelasticModel =
            new G4LEAntiProtonInelastic;
            theInelasticProcess->RegisterMe(theInelasticModel);
            theInelasticProcess->RegisterMe(theTheoModel);
            
            ChannelingDiscreteWrapper *theInelasticProcess_wrapper = new ChannelingDiscreteWrapper();
            theInelasticProcess_wrapper->RegisterProcess(theInelasticProcess);
            
            if(bWrapperOn){
                pManager->AddDiscreteProcess(theElasticProcess_wrapper);
                pManager->AddDiscreteProcess(theInelasticProcess_wrapper);
            }
            else{
                pManager->AddDiscreteProcess(theElasticProcess);
                pManager->AddDiscreteProcess(theInelasticProcess);
            }
        }
        else if (particleName == "neutron") {
            G4HadronElasticProcess* theElasticProcess1 =
            new G4HadronElasticProcess;

            // elastic scattering
            G4HadronElastic* theElasticModel1 = new G4HadronElastic;
            theElasticProcess1->RegisterMe(theElasticModel1);
            pManager->AddDiscreteProcess(theElasticProcess1);
            // inelastic scattering
            G4NeutronInelasticProcess* theInelasticProcess =
            new G4NeutronInelasticProcess("inelastic");
            G4LENeutronInelastic* theInelasticModel = new G4LENeutronInelastic;
            theInelasticProcess->RegisterMe(theInelasticModel);
            theInelasticProcess->RegisterMe(theTheoModel);
            
            pManager->AddDiscreteProcess(theElasticProcess1);
            pManager->AddDiscreteProcess(theInelasticProcess);
            
            // fission
            G4HadronFissionProcess* theFissionProcess =
            new G4HadronFissionProcess;
            G4LFission* theFissionModel = new G4LFission;
            theFissionProcess->RegisterMe(theFissionModel);
            pManager->AddDiscreteProcess(theFissionProcess);
            // capture
            G4HadronCaptureProcess* theCaptureProcess =
            new G4HadronCaptureProcess;
            G4NeutronRadCapture* theCaptureModel = new G4NeutronRadCapture;
            theCaptureProcess->RegisterMe(theCaptureModel);
            pManager->AddDiscreteProcess(theCaptureProcess);
        }
        else if (particleName == "anti_neutron") {
            
            G4AntiNeutronInelasticProcess* theInelasticProcess =
            new G4AntiNeutronInelasticProcess("inelastic");
            G4LEAntiNeutronInelastic* theInelasticModel =
            new G4LEAntiNeutronInelastic;
            theInelasticProcess->RegisterMe(theInelasticModel);
            theInelasticProcess->RegisterMe(theTheoModel);
            
            pManager->AddDiscreteProcess(theElasticProcess);
            pManager->AddDiscreteProcess(theInelasticProcess);
            
        }
        else if (particleName == "lambda") {
            
            G4LambdaInelasticProcess* theInelasticProcess =
            new G4LambdaInelasticProcess("inelastic");
            G4LELambdaInelastic* theInelasticModel = new G4LELambdaInelastic;
            theInelasticProcess->RegisterMe(theInelasticModel);
            theInelasticProcess->RegisterMe(theTheoModel);
            
            pManager->AddDiscreteProcess(theElasticProcess);
            pManager->AddDiscreteProcess(theInelasticProcess);
        }
        else if (particleName == "anti_lambda") {
            
            G4AntiLambdaInelasticProcess* theInelasticProcess =
            new G4AntiLambdaInelasticProcess("inelastic");
            G4LEAntiLambdaInelastic* theInelasticModel =
            new G4LEAntiLambdaInelastic;
            theInelasticProcess->RegisterMe(theInelasticModel);
            theInelasticProcess->RegisterMe(theTheoModel);
            
            pManager->AddDiscreteProcess(theElasticProcess);
            pManager->AddDiscreteProcess(theInelasticProcess);
            
        }
        else if (particleName == "sigma+") {
            
            G4SigmaPlusInelasticProcess* theInelasticProcess =
            new G4SigmaPlusInelasticProcess("inelastic");
            G4LESigmaPlusInelastic* theInelasticModel =
            new G4LESigmaPlusInelastic;
            theInelasticProcess->RegisterMe(theInelasticModel);
            theInelasticProcess->RegisterMe(theTheoModel);
            
            ChannelingDiscreteWrapper *theInelasticProcess_wrapper = new ChannelingDiscreteWrapper();
            theInelasticProcess_wrapper->RegisterProcess(theInelasticProcess);
            
            if(bWrapperOn){
                pManager->AddDiscreteProcess(theElasticProcess_wrapper);
                pManager->AddDiscreteProcess(theInelasticProcess_wrapper);
            }
            else{
                pManager->AddDiscreteProcess(theElasticProcess);
                pManager->AddDiscreteProcess(theInelasticProcess);
            }
        }
        else if (particleName == "sigma-") {
            
            G4SigmaMinusInelasticProcess* theInelasticProcess =
            new G4SigmaMinusInelasticProcess("inelastic");
            G4LESigmaMinusInelastic* theInelasticModel =
            new G4LESigmaMinusInelastic;
            theInelasticProcess->RegisterMe(theInelasticModel);
            theInelasticProcess->RegisterMe(theTheoModel);
            
            ChannelingDiscreteWrapper *theInelasticProcess_wrapper = new ChannelingDiscreteWrapper();
            theInelasticProcess_wrapper->RegisterProcess(theInelasticProcess);
            
            if(bWrapperOn){
                pManager->AddDiscreteProcess(theElasticProcess_wrapper);
                pManager->AddDiscreteProcess(theInelasticProcess_wrapper);
            }
            else{
                pManager->AddDiscreteProcess(theElasticProcess);
                pManager->AddDiscreteProcess(theInelasticProcess);
            }
        }
        else if (particleName == "anti_sigma+") {
            
            G4AntiSigmaPlusInelasticProcess* theInelasticProcess =
            new G4AntiSigmaPlusInelasticProcess("inelastic");
            G4LEAntiSigmaPlusInelastic* theInelasticModel =
            new G4LEAntiSigmaPlusInelastic;
            theInelasticProcess->RegisterMe(theInelasticModel);
            theInelasticProcess->RegisterMe(theTheoModel);
            
            ChannelingDiscreteWrapper *theInelasticProcess_wrapper = new ChannelingDiscreteWrapper();
            theInelasticProcess_wrapper->RegisterProcess(theInelasticProcess);
            
            if(bWrapperOn){
                pManager->AddDiscreteProcess(theElasticProcess_wrapper);
                pManager->AddDiscreteProcess(theInelasticProcess_wrapper);
            }
            else{
                pManager->AddDiscreteProcess(theElasticProcess);
                pManager->AddDiscreteProcess(theInelasticProcess);
            }
        }
        else if (particleName == "anti_sigma-") {
            
            G4AntiSigmaMinusInelasticProcess* theInelasticProcess =
            new G4AntiSigmaMinusInelasticProcess("inelastic");
            G4LEAntiSigmaMinusInelastic* theInelasticModel =
            new G4LEAntiSigmaMinusInelastic;
            theInelasticProcess->RegisterMe(theInelasticModel);
            theInelasticProcess->RegisterMe(theTheoModel);
            
            ChannelingDiscreteWrapper *theInelasticProcess_wrapper = new ChannelingDiscreteWrapper();
            theInelasticProcess_wrapper->RegisterProcess(theInelasticProcess);
            
            if(bWrapperOn){
                pManager->AddDiscreteProcess(theElasticProcess_wrapper);
                pManager->AddDiscreteProcess(theInelasticProcess_wrapper);
            }
            else{
                pManager->AddDiscreteProcess(theElasticProcess);
                pManager->AddDiscreteProcess(theInelasticProcess);
            }
        }
        else if (particleName == "xi0") {
            
            G4XiZeroInelasticProcess* theInelasticProcess =
            new G4XiZeroInelasticProcess("inelastic");
            G4LEXiZeroInelastic* theInelasticModel =
            new G4LEXiZeroInelastic;
            theInelasticProcess->RegisterMe(theInelasticModel);
            theInelasticProcess->RegisterMe(theTheoModel);
            
            pManager->AddDiscreteProcess(theElasticProcess);
            pManager->AddDiscreteProcess(theInelasticProcess);
            
        }
        else if (particleName == "xi-") {
            
            G4XiMinusInelasticProcess* theInelasticProcess =
            new G4XiMinusInelasticProcess("inelastic");
            G4LEXiMinusInelastic* theInelasticModel =
            new G4LEXiMinusInelastic;
            theInelasticProcess->RegisterMe(theInelasticModel);
            theInelasticProcess->RegisterMe(theTheoModel);
            
            ChannelingDiscreteWrapper *theInelasticProcess_wrapper = new ChannelingDiscreteWrapper();
            theInelasticProcess_wrapper->RegisterProcess(theInelasticProcess);
            
            if(bWrapperOn){
                pManager->AddDiscreteProcess(theElasticProcess_wrapper);
                pManager->AddDiscreteProcess(theInelasticProcess_wrapper);
            }
            else{
                pManager->AddDiscreteProcess(theElasticProcess);
                pManager->AddDiscreteProcess(theInelasticProcess);
            }
        }
        else if (particleName == "anti_xi0") {
            
            G4AntiXiZeroInelasticProcess* theInelasticProcess =
            new G4AntiXiZeroInelasticProcess("inelastic");
            G4LEAntiXiZeroInelastic* theInelasticModel =
            new G4LEAntiXiZeroInelastic;
            theInelasticProcess->RegisterMe(theInelasticModel);
            theInelasticProcess->RegisterMe(theTheoModel);
            
            pManager->AddDiscreteProcess(theElasticProcess);
            pManager->AddDiscreteProcess(theInelasticProcess);
            
        }
        else if (particleName == "anti_xi-") {
            
            G4AntiXiMinusInelasticProcess* theInelasticProcess =
            new G4AntiXiMinusInelasticProcess("inelastic");
            G4LEAntiXiMinusInelastic* theInelasticModel =
            new G4LEAntiXiMinusInelastic;
            theInelasticProcess->RegisterMe(theInelasticModel);
            theInelasticProcess->RegisterMe(theTheoModel);
            
            ChannelingDiscreteWrapper *theInelasticProcess_wrapper = new ChannelingDiscreteWrapper();
            theInelasticProcess_wrapper->RegisterProcess(theInelasticProcess);
            
            if(bWrapperOn){
                pManager->AddDiscreteProcess(theElasticProcess_wrapper);
                pManager->AddDiscreteProcess(theInelasticProcess_wrapper);
            }
            else{
                pManager->AddDiscreteProcess(theElasticProcess);
                pManager->AddDiscreteProcess(theInelasticProcess);
            }
        }
        else if (particleName == "deuteron") {
            
            G4DeuteronInelasticProcess* theInelasticProcess =
            new G4DeuteronInelasticProcess("inelastic");
            G4LEDeuteronInelastic* theInelasticModel =
            new G4LEDeuteronInelastic;
            theInelasticProcess->RegisterMe(theInelasticModel);
            theInelasticProcess->RegisterMe(theTheoModel);
            
            pManager->AddDiscreteProcess(theElasticProcess);
            pManager->AddDiscreteProcess(theInelasticProcess);
            
        }
        else if (particleName == "triton") {
            
            G4TritonInelasticProcess* theInelasticProcess =
            new G4TritonInelasticProcess("inelastic");
            G4LETritonInelastic* theInelasticModel =
            new G4LETritonInelastic;
            theInelasticProcess->RegisterMe(theInelasticModel);
            theInelasticProcess->RegisterMe(theTheoModel);
            
            pManager->AddDiscreteProcess(theElasticProcess);
            pManager->AddDiscreteProcess(theInelasticProcess);
            
        }
        else if (particleName == "alpha") {
            
            G4AlphaInelasticProcess* theInelasticProcess =
            new G4AlphaInelasticProcess("inelastic");
            G4LEAlphaInelastic* theInelasticModel =
            new G4LEAlphaInelastic;
            theInelasticProcess->RegisterMe(theInelasticModel);
            theInelasticProcess->RegisterMe(theTheoModel);
            
            ChannelingDiscreteWrapper *theInelasticProcess_wrapper = new ChannelingDiscreteWrapper();
            theInelasticProcess_wrapper->RegisterProcess(theInelasticProcess);
            
            if(bWrapperOn){
                pManager->AddDiscreteProcess(theElasticProcess_wrapper);
                pManager->AddDiscreteProcess(theInelasticProcess_wrapper);
            }
            else{
                pManager->AddDiscreteProcess(theElasticProcess);
                pManager->AddDiscreteProcess(theInelasticProcess);
            }
        }
        else if (particleName == "omega-") {
            
            G4OmegaMinusInelasticProcess* theInelasticProcess =
            new G4OmegaMinusInelasticProcess("inelastic");
            G4LEOmegaMinusInelastic* theInelasticModel =
            new G4LEOmegaMinusInelastic;
            theInelasticProcess->RegisterMe(theInelasticModel);
            theInelasticProcess->RegisterMe(theTheoModel);
            
            ChannelingDiscreteWrapper *theInelasticProcess_wrapper = new ChannelingDiscreteWrapper();
            theInelasticProcess_wrapper->RegisterProcess(theInelasticProcess);
            
            if(bWrapperOn){
                pManager->AddDiscreteProcess(theElasticProcess_wrapper);
                pManager->AddDiscreteProcess(theInelasticProcess_wrapper);
            }
            else{
                pManager->AddDiscreteProcess(theElasticProcess);
                pManager->AddDiscreteProcess(theInelasticProcess);
            }
        }
        else if (particleName == "anti_omega-") {
            
            G4AntiOmegaMinusInelasticProcess* theInelasticProcess =
            new G4AntiOmegaMinusInelasticProcess("inelastic");
            G4LEAntiOmegaMinusInelastic* theInelasticModel =
            new G4LEAntiOmegaMinusInelastic;
            theInelasticProcess->RegisterMe(theInelasticModel);
            theInelasticProcess->RegisterMe(theTheoModel);
            
            ChannelingDiscreteWrapper *theInelasticProcess_wrapper = new ChannelingDiscreteWrapper();
            theInelasticProcess_wrapper->RegisterProcess(theInelasticProcess);
            
            if(bWrapperOn){
                pManager->AddDiscreteProcess(theElasticProcess_wrapper);
                pManager->AddDiscreteProcess(theInelasticProcess_wrapper);
            }
            else{
                pManager->AddDiscreteProcess(theElasticProcess);
                pManager->AddDiscreteProcess(theInelasticProcess);
            }
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChannelingPhysicsList::AddStandardSS(G4ParticleDefinition* particle){
    G4ProcessManager* pManager = particle->GetProcessManager();
    
    G4CoulombScattering* cs = new G4CoulombScattering();
    cs->SetBuildTableFlag(false);
    ChannelingDiscreteWrapper *cs_wrapper = new ChannelingDiscreteWrapper();
    cs_wrapper->RegisterProcess(cs,1);
    
    if(bWrapperOn){
        pManager->AddDiscreteProcess(cs_wrapper,1);
    }
    else{
        pManager->AddDiscreteProcess(cs,1);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChannelingPhysicsList::AddStandardNR(G4ParticleDefinition* particle){
    G4ProcessManager* pManager = particle->GetProcessManager();
    
    G4ScreenedNuclearRecoil* nucr = new G4ScreenedNuclearRecoil();
    ChannelingDiscreteWrapper *nucr_wrapper = new ChannelingDiscreteWrapper();
    nucr_wrapper->RegisterProcess(nucr,1);
    
    if(bWrapperOn){
        pManager->AddDiscreteProcess(nucr_wrapper,1);
    }
    else{
        pManager->AddDiscreteProcess(nucr,1);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChannelingPhysicsList::AddEM()
{    
    theParticleIterator->reset();
    while( (*theParticleIterator)() ){
        G4ParticleDefinition* particle = theParticleIterator->value();
        G4ProcessManager* pManager = particle->GetProcessManager();
        G4String particleName = particle->GetParticleName();
        
        if (particleName == "gamma") {
            // gamma
            // Construct processes for gamma
            // Compton scattering
            G4ComptonScattering* cs = new G4ComptonScattering;
            cs->SetEmModel(new G4KleinNishinaModel(),1);
            G4VEmModel* theLowEPComptonModel = new G4LowEPComptonModel();
            theLowEPComptonModel->SetHighEnergyLimit(20*MeV);
            cs->AddEmModel(0, theLowEPComptonModel);
            
            // Photoelectric
            G4PhotoElectricEffect* pe = new G4PhotoElectricEffect();
            G4VEmModel* theLivermorePEModel = new G4LivermorePhotoElectricModel();
            theLivermorePEModel->SetHighEnergyLimit(10*GeV);
            pe->SetEmModel(theLivermorePEModel,1);
            
            // Gamma conversion
            G4GammaConversion* gc = new G4GammaConversion();
            G4VEmModel* thePenelopeGCModel = new G4PenelopeGammaConversionModel();
            thePenelopeGCModel->SetHighEnergyLimit(1*GeV);
            gc->SetEmModel(thePenelopeGCModel,1);
            
            // Rayleigh scattering
            G4RayleighScattering* rs = new G4RayleighScattering();
            
            pManager->AddDiscreteProcess(cs);
            pManager->AddDiscreteProcess(pe);
            pManager->AddDiscreteProcess(gc);
            pManager->AddDiscreteProcess(rs);
            
        }
        else if (particleName == "e-") {
            //electron
            // Construct processes for electron
            // single scattering
            G4CoulombScattering* ecs = new G4CoulombScattering();
            ecs->SetBuildTableFlag(false);
            ChannelingDiscreteWrapper *ecs_wrapper = new ChannelingDiscreteWrapper();
            ecs_wrapper->RegisterProcess(ecs,1);
            
            // ionisation
            G4eIonisation* eIoni = new G4eIonisation();
            eIoni->SetStepFunction(0.2, 100*um);
            G4VEmModel* theIoniPenelope = new G4PenelopeIonisationModel();
            theIoniPenelope->SetHighEnergyLimit(0.1*MeV);
            eIoni->AddEmModel(0, theIoniPenelope, new G4UniversalFluctuation());
            
            ChannelingContinuousDiscreteWrapper *eIoni_wrapper = new ChannelingContinuousDiscreteWrapper();
            eIoni_wrapper->RegisterProcess(eIoni,-1);
            
            // bremsstrahlung
            G4eBremsstrahlung* eBrem = new G4eBremsstrahlung();
            
            ChannelingContinuousDiscreteWrapper *eBrem_wrapper = new ChannelingContinuousDiscreteWrapper();
            eBrem_wrapper->RegisterProcess(eBrem,-1);
            
            if(bWrapperOn){
                pManager->AddDiscreteProcess(ecs_wrapper,1);
                pManager->AddProcess(eIoni_wrapper,-1, 2, 2);
                pManager->AddProcess(eBrem_wrapper,-1,-1, 3);
                
            }
            else{
                pManager->AddDiscreteProcess(ecs,1);
                pManager->AddProcess(eIoni,-1, 2, 2);
                pManager->AddProcess(eBrem,-1,-1, 3);
            }
        }
        else if (particleName == "e+") {
            //positron
            // Construct processes for positron
            // single scattering
            G4CoulombScattering* ecs = new G4CoulombScattering();
            ecs->SetBuildTableFlag(false);
            ChannelingDiscreteWrapper *ecs_wrapper = new ChannelingDiscreteWrapper();
            ecs_wrapper->RegisterProcess(ecs,1);
            
            // ionisation
            G4eIonisation* eIoni = new G4eIonisation();
            eIoni->SetStepFunction(0.2, 100*um);
            G4VEmModel* theIoniPenelope = new G4PenelopeIonisationModel();
            theIoniPenelope->SetHighEnergyLimit(0.1*MeV);
            eIoni->AddEmModel(0, theIoniPenelope, new G4UniversalFluctuation());
            
            ChannelingContinuousDiscreteWrapper *eIoni_wrapper = new ChannelingContinuousDiscreteWrapper();
            eIoni_wrapper->RegisterProcess(eIoni,-1);
            
            // bremsstrahlung
            G4eBremsstrahlung* eBrem = new G4eBremsstrahlung();
            
            ChannelingContinuousDiscreteWrapper *eBrem_wrapper = new ChannelingContinuousDiscreteWrapper();
            eBrem_wrapper->RegisterProcess(eBrem,-1);
            
            // annihilation
            G4eplusAnnihilation* eplusAnn = new G4eplusAnnihilation();
            
            ChannelingDiscreteWrapper *eplusAnn_wrapper = new ChannelingDiscreteWrapper();
            eplusAnn_wrapper->RegisterProcess(eplusAnn,-1);
            
            if(bWrapperOn){
                pManager->AddDiscreteProcess(ecs_wrapper,1);
                pManager->AddProcess(eIoni_wrapper,-1, 2, 2);
                pManager->AddProcess(eBrem_wrapper,-1, -1, 3);
                //pManager->AddProcess(eplusAnn_wrapper, 0,-1, 4);
                pManager->AddDiscreteProcess(eplusAnn_wrapper, 0);
            }
            else{
                pManager->AddDiscreteProcess(ecs,1);
                pManager->AddProcess(eIoni,-1, 2, 2);
                pManager->AddProcess(eBrem,-1,-1, 3);
                pManager->AddProcess(eplusAnn, 0,-1, 4);
            }
        }
        else if( particleName == "mu+" ||
                particleName == "mu-"    ) {
            //muon
            // Construct processes for muon+
            if(fScatteringType == "ss"){
                AddStandardSS(particle);
            }
            else if(fScatteringType == "nr"){
                AddStandardNR(particle);
            }
            
            // ionisation
            G4MuIonisation* muIoni = new G4MuIonisation();
            
            ChannelingContinuousDiscreteWrapper *muIoni_wrapper = new ChannelingContinuousDiscreteWrapper();
            muIoni_wrapper->RegisterProcess(muIoni,-1);
            
            // bremsstrahlung
            G4MuBremsstrahlung* muBrem = new G4MuBremsstrahlung();
            
            ChannelingContinuousDiscreteWrapper *muBrem_wrapper = new ChannelingContinuousDiscreteWrapper();
            muBrem_wrapper->RegisterProcess(muBrem,-1);
            
            // annihilation
            G4MuPairProduction* muPairProd = new G4MuPairProduction();
            
            ChannelingContinuousDiscreteWrapper* muPairProd_wrapper = new ChannelingContinuousDiscreteWrapper();
            muPairProd_wrapper->RegisterProcess(muPairProd,-1);
            
            if(bWrapperOn){
                pManager->AddProcess(muIoni_wrapper,-1, 2, 2);
                pManager->AddProcess(muBrem_wrapper,-1, -1, 3);
                pManager->AddProcess(muPairProd_wrapper,-1,-1, 4);
            }
            else{
                pManager->AddProcess(muIoni,-1, 2, 2);
                pManager->AddProcess(muBrem,-1,-1, 3);
                pManager->AddProcess(muPairProd,-1,-1, 4);
            }
        }
        else if( particleName == "GenericIon" ) {
            if(fScatteringType == "ss"){
                AddStandardSS(particle);
            }
            else if(fScatteringType == "nr"){
                AddStandardNR(particle);
            }
            
            // nuclear stopping
            G4NuclearStopping* ionnuc = new G4NuclearStopping();
            ChannelingDiscreteWrapper *ionnuc_wrapper = new ChannelingDiscreteWrapper();
            ionnuc->SetMaxKinEnergy(MeV);
            ionnuc_wrapper->RegisterProcess(ionnuc,+1);

            // ionisation
            G4hIonisation* theppIonisation = new G4hIonisation();
            ChannelingContinuousDiscreteWrapper *theppIonisation_wrapper = new ChannelingContinuousDiscreteWrapper();
            theppIonisation_wrapper->RegisterProcess(theppIonisation,-1);
            
            if(bWrapperOn){
                pManager->AddProcess(theppIonisation_wrapper,-1, 2, 2);
                pManager->AddProcess(ionnuc_wrapper);
            }
            else{
                pManager->AddProcess(theppIonisation,-1, 2, 2);
                pManager->AddProcess(ionnuc);
            }
        }
        else {
            if ((particle->GetPDGCharge() != 0.0) &&
                (particle->GetParticleName() != "chargedgeantino")&&
                (!particle->IsShortLived()) ) {
                if(fScatteringType == "ss"){
                    AddStandardSS(particle);
                }
                else if(fScatteringType == "nr"){
                    AddStandardNR(particle);
                }
                
                // nuclear stopping
                G4NuclearStopping* ionnuc = new G4NuclearStopping();
                ChannelingDiscreteWrapper *ionnuc_wrapper = new ChannelingDiscreteWrapper();
                ionnuc->SetMaxKinEnergy(MeV);
                ionnuc_wrapper->RegisterProcess(ionnuc,+1);

                // ionisation
                G4hIonisation* theppIonisation = new G4hIonisation();
                ChannelingContinuousDiscreteWrapper *theppIonisation_wrapper = new ChannelingContinuousDiscreteWrapper();
                theppIonisation_wrapper->RegisterProcess(theppIonisation,-1);
                
                if(bWrapperOn){
                    pManager->AddProcess(theppIonisation_wrapper,-1, 2, 2);
                    pManager->AddProcess(ionnuc_wrapper);
                }
                else{
                    pManager->AddProcess(theppIonisation,-1, 2, 2);
                    pManager->AddProcess(ionnuc);
               }
            }
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChannelingPhysicsList::AddChanneling(){
    
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
    
    XCrystalIntegratedDensityHub* vIntegratedDensityHub = new XCrystalIntegratedDensityHub();
    vIntegratedDensityHub->SetPotential(vPotentialEnergy);
    vIntegratedDensityHub->SetDensityNuclei(vNucleiDensity);
    vIntegratedDensityHub->SetDensityElectron(vElectronDensity);
    
    for(G4int i=-3;i<=+3;i++){
        if(i!=0){
            vIntegratedDensityHub->SetIntegratedDensityNuclei(new XCrystalIntegratedDensityPlanar(),i);
            vIntegratedDensityHub->SetIntegratedDensityElectron(new XCrystalIntegratedDensityPlanar(),i);
        }
    }
    
    G4CMPChanneling* channeling =  new G4CMPChanneling();
    //G4CMPChannelingSlow* channeling =  new G4CMPChannelingSlow();
    channeling->SetPotential(vPotentialEnergy);
    channeling->SetIntegratedDensity(vIntegratedDensityHub);
    channeling->SetElectricField(vElectricField);
    channeling->SetFileName(fFileName);
    
    theParticleIterator->reset();
    while( (*theParticleIterator)() ){
        G4ParticleDefinition* particle = theParticleIterator->value();
        G4ProcessManager* pManager = particle->GetProcessManager();
        G4String particleName = particle->GetParticleName();
        
        if(particle->GetPDGCharge() != 0)
        {
            pManager->AddDiscreteProcess(channeling);
        }
    }
    G4cout<<"\nPhysicsList::AddChanneling: Channeling process added...\n"<<G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChannelingPhysicsList::AddDecay(){
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
    G4cout<<"\nPhysicsList::AddDecay: Decay processes added...\n"<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChannelingPhysicsList::AddChannelingStepLimiter()
{
    // Step limitation seen as a process
    ChannelingStepLimiter* stepMaxProcess = new ChannelingStepLimiter();
    
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


void ChannelingPhysicsList::SetCuts()
{
    // These values are used as the default production thresholds
    // for the world volume.
    SetCutsWithDefault();
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ChannelingPhysicsList::SetFileName(const G4String& vFilename){
    if(fFileName != vFilename){
        fFileName = vFilename;
        G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4String ChannelingPhysicsList::GetFileName(){
    return fFileName;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ChannelingPhysicsList::SetScatteringType(const G4String& vScatteringType){
    if(fScatteringType != vScatteringType){
        fScatteringType = vScatteringType;
        G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4String ChannelingPhysicsList::GetScatteringType(){
    return fScatteringType;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ChannelingPhysicsList::EnableChanneling(G4bool flag) {
    if(bChannelingOn != flag){
        bChannelingOn = flag;
        G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ChannelingPhysicsList::EnableWrapper(G4bool flag) {
    if(bWrapperOn != flag){
        bWrapperOn = flag;
        G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ChannelingPhysicsList::EnableDecay(G4bool flag) {
    if(bDecayOn != flag){
        bDecayOn = flag;
        G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


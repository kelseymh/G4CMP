/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#include "ChannelingContinuousDiscreteWrapper.hh"

#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4VParticleChange.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "G4RandomTools.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "G4GeometryTolerance.hh"

#include "G4SystemOfUnits.hh"


ChannelingContinuousDiscreteWrapper::ChannelingContinuousDiscreteWrapper(const G4String& aName)
:G4VContinuousDiscreteProcess(aName){
    if (verboseLevel>1) {
        G4cout << GetProcessName() << " is created "<< G4endl;
    }
    bNucleiOrElectronFlag = +0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ChannelingContinuousDiscreteWrapper::ChannelingContinuousDiscreteWrapper(const G4String&, G4VContinuousDiscreteProcess* toRegister)
:G4VContinuousDiscreteProcess("ChannelingContinuousDiscreteWrapper"){
    fRegisteredProcess = toRegister;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ChannelingContinuousDiscreteWrapper::~ChannelingContinuousDiscreteWrapper(){
    ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ChannelingContinuousDiscreteWrapper::ChannelingContinuousDiscreteWrapper(ChannelingContinuousDiscreteWrapper& right): G4VContinuousDiscreteProcess(right){
    ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ChannelingContinuousDiscreteWrapper::RegisterProcess(G4VContinuousDiscreteProcess* toRegister){
    fRegisteredProcess = toRegister;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ChannelingContinuousDiscreteWrapper::RegisterProcess(G4VContinuousDiscreteProcess* toRegister, G4int flag){
    fRegisteredProcess = toRegister;
    bNucleiOrElectronFlag = flag;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ChannelingContinuousDiscreteWrapper::SetNucleiOrElectronFlag(G4int flag){
    bNucleiOrElectronFlag = flag;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int ChannelingContinuousDiscreteWrapper::GetNucleiOrElectronFlag(){
    return bNucleiOrElectronFlag;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double ChannelingContinuousDiscreteWrapper::GetDensity(const G4Track& aTrack){
    //Retrieve nuclei and electron density from ChannelingParticleUserInfo object
    ChannelingParticleUserInfo* chanInfo = (ChannelingParticleUserInfo*) aTrack.GetUserInformation();
    
    G4double vDensity = 1.;
    
    if(chanInfo){
        if(bNucleiOrElectronFlag == +1){
            vDensity = chanInfo->GetNucleiDensity();
        }
        else if(bNucleiOrElectronFlag == -1){
            vDensity = chanInfo->GetElectronDensity();
        }
        else{
            vDensity = (chanInfo->GetNucleiDensity() + chanInfo->GetElectronDensity())/2.;
        }
    }
    else {
        G4cout << G4endl << "ChannelingContinuousDiscreteWrapper:: ERROR - no ChannelingParticleUserInfo object Detected" << G4endl;
    }
    
    return vDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double ChannelingContinuousDiscreteWrapper::GetDensityPreviousStep(const G4Track& aTrack){
    //Retrieve nuclei and electron density from ChannelingParticleUserInfo object
    ChannelingParticleUserInfo* chanInfo = (ChannelingParticleUserInfo*) aTrack.GetUserInformation();
    
    G4double vDensityPreviousStep = 1.;
    
    if(chanInfo){
        if(bNucleiOrElectronFlag == +1){
            vDensityPreviousStep = chanInfo->GetNucleiDensityPreviousStep();
        }
        else if(bNucleiOrElectronFlag == -1){
            vDensityPreviousStep = chanInfo->GetElectronDensityPreviousStep();
        }
        else{
            vDensityPreviousStep = (chanInfo->GetNucleiDensityPreviousStep() + chanInfo->GetElectronDensityPreviousStep())/2.;
        }
    }
    else {
        G4cout << G4endl << "ChannelingContinuousDiscreteWrapper:: ERROR - no ChannelingParticleUserInfo object Detected" << G4endl;
    }
    
    return vDensityPreviousStep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ChannelingContinuousDiscreteWrapper::StartTracking(G4Track* aTrack){
    fRegisteredProcess->StartTracking(aTrack);
    currentInteractionLength = -1.0;
    theNumberOfInteractionLengthLeft = -1.0;
    theInitialNumberOfInteractionLength = -1.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double ChannelingContinuousDiscreteWrapper::GetMeanFreePath(const G4Track&,
                                                            G4double, //previousStepSize,
                                                            G4ForceCondition*){
    return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double ChannelingContinuousDiscreteWrapper::PostStepGetPhysicalInteractionLength (const G4Track &aTrack,
                                                                                  G4double previousStepSize,
                                                                                  G4ForceCondition *condition){
    
    G4double vDensity = GetDensity(aTrack);
    G4double vDensityPreviousStep = GetDensityPreviousStep(aTrack);
    
    if ( (previousStepSize < 0.0) || (theNumberOfInteractionLengthLeft<=0.0)) {
        // beginning of tracking (or just after DoIt of this process)
        ResetNumberOfInteractionLengthLeft();
    } else if ( previousStepSize > 0.0) {
        // subtract NumberOfInteractionLengthLeft
        SubtractNumberOfInteractionLengthLeft(previousStepSize * vDensityPreviousStep);
    } else {
        // zero step DO NOTHING
    }
    
    G4double regIntLength = fRegisteredProcess->PostStepGetPhysicalInteractionLength(aTrack, previousStepSize * vDensityPreviousStep, condition);
    G4double regIntNumber = fRegisteredProcess->GetNumberOfInteractionLengthLeft();
    currentInteractionLength = regIntLength / regIntNumber;
    theNumberOfInteractionLengthLeft = regIntNumber;
    
    currentInteractionLength =  theNumberOfInteractionLengthLeft * currentInteractionLength;
    currentInteractionLength /= vDensity;
    return currentInteractionLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double ChannelingContinuousDiscreteWrapper::AlongStepGetPhysicalInteractionLength (const G4Track& aTrack,
                                                                                   G4double  previousStepSize,
                                                                                   G4double  currentMinimumStep,
                                                                                   G4double& currentSafety,
                                                                                   G4GPILSelection* selection){
    G4double vDensityPreviousStep = GetDensityPreviousStep(aTrack);

    return fRegisteredProcess->AlongStepGetPhysicalInteractionLength(aTrack,previousStepSize * vDensityPreviousStep,currentMinimumStep,currentSafety,selection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* ChannelingContinuousDiscreteWrapper::PostStepDoIt(const G4Track& aTrack,
                                                                   const G4Step& aStep ){
    G4double vDensity = GetDensity(aTrack);
    G4double vStepLengthSaved = aStep.GetStepLength();
    const_cast<G4Step&>(aStep).SetStepLength(aStep.GetStepLength() * vDensity);
    G4VParticleChange* theParticleChange = fRegisteredProcess->PostStepDoIt(aTrack, aStep);
    const_cast<G4Step&>(aStep).SetStepLength(vStepLengthSaved);
    return theParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* ChannelingContinuousDiscreteWrapper::AlongStepDoIt(const G4Track& aTrack,
                                                                    const G4Step& aStep ){
    G4double vDensity = GetDensity(aTrack);
    G4double vStepLengthSaved = aStep.GetStepLength();
    const_cast<G4Step&>(aStep).SetStepLength(aStep.GetStepLength() * vDensity);
    G4VParticleChange* theParticleChange = fRegisteredProcess->AlongStepDoIt(aTrack, aStep);
    const_cast<G4Step&>(aStep).SetStepLength(vStepLengthSaved);
    return theParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4bool ChannelingContinuousDiscreteWrapper::IsApplicable(const G4ParticleDefinition& aParticleDefinition){
    return fRegisteredProcess->IsApplicable(aParticleDefinition);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ChannelingContinuousDiscreteWrapper::BuildPhysicsTable(const G4ParticleDefinition& aParticleDefinition){
    fRegisteredProcess->BuildPhysicsTable(aParticleDefinition);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ChannelingContinuousDiscreteWrapper::PreparePhysicsTable(const G4ParticleDefinition& aParticleDefinition){
    fRegisteredProcess->PreparePhysicsTable(aParticleDefinition);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool ChannelingContinuousDiscreteWrapper::StorePhysicsTable(const G4ParticleDefinition* aParticleDefinition,
                                                            const G4String& aString,
                                                            G4bool aBool){
    return fRegisteredProcess->StorePhysicsTable(aParticleDefinition,aString,aBool);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool ChannelingContinuousDiscreteWrapper::RetrievePhysicsTable( const G4ParticleDefinition* aParticleDefinition,
                                                               const G4String& aString,
                                                               G4bool aBool){
    return fRegisteredProcess->RetrievePhysicsTable(aParticleDefinition,aString,aBool);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double ChannelingContinuousDiscreteWrapper::GetContinuousStepLimit(const G4Track& ,G4double  ,G4double ,G4double& ){
    return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

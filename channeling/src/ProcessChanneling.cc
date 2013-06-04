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

#include "ProcessChanneling.hh"

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

#include "XLatticeManager3.hh"

#include "XLogicalAtomicLattice.hh"
#include "XLogicalAtomicLatticeDiamond.hh"
#include "XLogicalBase.hh"
#include "XUnitCell.hh"

#include "XCrystalPlanarMoliereTempPotential.hh"
#include "XCrystalPlanarMoliereElectricField.hh"
#include "XCrystalPlanarNucleiDensity.hh"
#include "XCrystalPlanarMoliereElectronDensity.hh"

#include "XCrystalIntegratedDensity.hh"

ProcessChanneling::ProcessChanneling(const G4String& aName):G4VDiscreteProcess(aName){
    fLatticeManager = XLatticeManager3::GetXLatticeManager();
    
    G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
    G4cout<<"\n ProcessChanneling::Constructor: Geometry surface tolerance is: " << kCarTolerance / mm << " mm";
    if(verboseLevel>1) {
        G4cout << GetProcessName() << " is created "<< G4endl;
    }
    
    fCompute = true;
    fHasBeenInChanneling = false;
    fFileOut.open("channeling.txt");

    InitializeCrystalCharacteristics();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ProcessChanneling::InitializeCrystalCharacteristics(){
    fPotentialEnergy = new XCrystalPlanarMoliereTempPotential();
    fElectricField = new XCrystalPlanarMoliereElectricField();
    fNucleiDensity = new XCrystalPlanarNucleiDensity();
    fElectronDensity = new XCrystalPlanarMoliereElectronDensity();
    fIntegratedDensity = new XCrystalIntegratedDensity(fPotentialEnergy,fNucleiDensity,fElectronDensity);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ProcessChanneling::~ProcessChanneling(){
    fFileOut.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ProcessChanneling::ProcessChanneling(ProcessChanneling& right):G4VDiscreteProcess(right){
    ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XVCrystalCharacteristic* ProcessChanneling::GetPotential(){
    return fPotentialEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ProcessChanneling::SetPotential(XVCrystalCharacteristic* vPotential){
    fPotentialEnergy = vPotential;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XVCrystalCharacteristic* ProcessChanneling::GetNucleiDensity(){
    return fNucleiDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ProcessChanneling::SetNucleiDensity(XVCrystalCharacteristic* vNucleiDensity){
    fNucleiDensity = vNucleiDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XVCrystalCharacteristic* ProcessChanneling::GetElectronDensity(){
    return fElectronDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ProcessChanneling::SetElectronDensity(XVCrystalCharacteristic* vElectronDensity){
    fElectronDensity = vElectronDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XVCrystalCharacteristic* ProcessChanneling::GetElectricField(){
    return fElectricField;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ProcessChanneling::SetElectricField(XVCrystalCharacteristic* vElectricField){
    fElectricField = vElectricField;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XCrystalIntegratedDensity* ProcessChanneling::GetIntegratedDensity(){
    return fIntegratedDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ProcessChanneling::SetIntegratedDensity(XCrystalIntegratedDensity* vIntegratedDensity){
    fIntegratedDensity = vIntegratedDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ProcessChanneling::UpdatePositionMomentumDensity(const G4Track& aTrack){
    XPhysicalLattice* vXtalPhysLattice = fLatticeManager->GetXPhysicalLattice(aTrack.GetStep()->GetPostStepPoint()->GetPhysicalVolume());
    G4VPhysicalVolume* vVolume = aTrack.GetStep()->GetPostStepPoint()->GetPhysicalVolume();
    
    if(!fIntegratedDensity->HasBeenInitialized()){
        fIntegratedDensity->SetXPhysicalLattice(fLatticeManager->GetXPhysicalLattice(vVolume));
        fIntegratedDensity->InitializeTable();
        G4cout << "ChannelingProcess::UpdatePositionMomentumDensity::fIntegratedDensity->Initialized" << std::endl;
    }

    switch (fHasBeenInChanneling) {
        case 0:
            fPosition.setX(G4UniformRand() * vXtalPhysLattice->ComputeInterplanarPeriod());
            fMomentum = fLatticeManager->GetXPhysicalLattice(vVolume)->ProjectVectorFromWorldToLattice(aTrack.GetMomentum());
            fDensity = fIntegratedDensity->GetValue(ComputeTransverseEnergy(aTrack).x());
            break;
        default:
            fMomentum += fLatticeManager->GetXPhysicalLattice(vVolume)->ProjectVectorFromWorldToLattice(aTrack.GetMomentum());
            fDensity = fIntegratedDensity->GetValue(ComputeTransverseEnergy(aTrack).x());
            break;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool ProcessChanneling::IsInChanneling(const G4Track& aTrack){
    //----------------------------------------
    // check if the particle momentum
    // transverse to the (h,k,l) plane
    // is small enough to permit channeling
    //----------------------------------------
    
    UpdatePositionMomentumDensity(aTrack);
    if( ComputeTransverseEnergy(aTrack).x() < ComputeChannelingCriticalEnergy(aTrack) ){
        fHasBeenInChanneling = true;
        return true;
    }
    return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double ProcessChanneling::ComputeChannelingCriticalEnergy(const G4Track& aTrack){
    
    XPhysicalLattice* vXtalPhysLattice = fLatticeManager->GetXPhysicalLattice(aTrack.GetStep()->GetPostStepPoint()->GetPhysicalVolume());
    
    G4double vCriticalEnergy = 0.;
    vCriticalEnergy += fPotentialEnergy->GetMaximum(vXtalPhysLattice).x();
    
    vCriticalEnergy -= fPotentialEnergy->GetMinimum(vXtalPhysLattice).x();
    
    return vCriticalEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector ProcessChanneling::ComputeTransverseEnergy(const G4Track& aTrack){
    
    XPhysicalLattice* vXtalPhysLattice = fLatticeManager->GetXPhysicalLattice(aTrack.GetStep()->GetPostStepPoint()->GetPhysicalVolume());
    
    G4ThreeVector vTransverseEnergy = G4ThreeVector(0.,0.,0.);
    
    G4double vTotalEnergy = aTrack.GetStep()->GetPostStepPoint()->GetTotalEnergy();
    
    vTransverseEnergy += fPotentialEnergy->ComputeValue(fPosition,vXtalPhysLattice);
        
    vTransverseEnergy += G4ThreeVector( 0., 0., (  pow( fMomentum.z(), 2. )  / vTotalEnergy ) );
    
    vTransverseEnergy += G4ThreeVector( (  pow( fMomentum.x(), 2. ) / vTotalEnergy ), 0., 0. );

    return vTransverseEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double
ProcessChanneling::GetChannelingMeanFreePath(const G4Track& aTrack)
{
//    // dechanneling length is the dacay length for the dechanneling processes
//    G4double vChannelingMeanFreePathNearNuclei = 1.5 * mm; // dechanneling length for particles which enter the crystal near nuclei
//    G4double vChannelingMeanFreePathFarFromNuclei = 20. * cm; // dechannelign length for particles which enter the crystal far from nuclei
//    G4double vParticleFractionNearNuclei = 0.2; // fraction of particles which enter the crystal near the nuclei
//    XPhysicalLattice* vXtalPhysLattice = fLatticeManager->GetXPhysicalLattice(aTrack.GetStep()->GetPostStepPoint()->GetPhysicalVolume());
//    
//    if(fPosition.x() / vXtalPhysLattice->ComputeInterplanarPeriod() < vParticleFractionNearNuclei){
//        return vChannelingMeanFreePathNearNuclei;
//    }
//    else{
//        return vChannelingMeanFreePathFarFromNuclei;
//    }
    return 60.E-6 * meter;
//    return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double
ProcessChanneling::GetMeanFreePath(const G4Track& aTrack, G4double /*previousStepSize*/, G4ForceCondition* condition  )
{
    //----------------------------------------
    // the condition is forced to check if
    // the volume has a lattice at each step
    //----------------------------------------
    
    *condition = Forced;
    
    G4VPhysicalVolume* vVolume = aTrack.GetStep()->GetPostStepPoint()->GetPhysicalVolume();
    
    if(fLatticeManager->HasLattice(vVolume))
    {
        if(IsInChanneling(aTrack)){
            return GetChannelingMeanFreePath(aTrack);
        }
    }
    else{
        return DBL_MAX;
    }
    
    return DBL_MAX;
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* ProcessChanneling::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep){
    
    aParticleChange.Initialize(aTrack);
    G4VPhysicalVolume* vVolume = aTrack.GetStep()->GetPostStepPoint()->GetPhysicalVolume();
    bool bIsInChanneling = false;
    
    if(fLatticeManager->HasLattice(vVolume))
    {
        if(bIsInChanneling = IsInChanneling(aTrack)){
            aParticleChange.ProposeMomentumDirection(fLatticeManager->GetXPhysicalLattice(vVolume)->GetLatticeDirection().unit());
        }
    }
    
    if(!bIsInChanneling)
    {
        if(fHasBeenInChanneling){
            aParticleChange.ProposeMomentumDirection(fLatticeManager->GetXPhysicalLattice(vVolume)->GetLatticeDirection().unit());
            fFileOut << bIsInChanneling << " " << fHasBeenInChanneling << " " << aTrack.GetStep()->GetPostStepPoint()->GetPosition().y() << " " << fPosition.x() << " " << fMomentum.x() << " " << fDensity << " " <<  ComputeTransverseEnergy(aTrack) <<  std::endl;
        }
        fHasBeenInChanneling = false;
    }
    

    return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4bool ProcessChanneling::IsApplicable(const G4ParticleDefinition& aPD)
{
    return(aPD.GetPDGCharge()>0.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void ProcessChanneling::BuildPhysicsTable(const G4ParticleDefinition&)
{
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void ProcessChanneling::ComputeCrystalCharacteristicForChanneling(const G4Track& aTrack){
    
    XPhysicalLattice* vXtalPhysLattice = fLatticeManager->GetXPhysicalLattice(aTrack.GetStep()->GetPostStepPoint()->GetPhysicalVolume());
    
    std::ofstream vFileOutPot;
    std::ofstream vFileOutEfx;
    std::ofstream vFileOutNud;
    std::ofstream vFileOutEld;
    
    vFileOutPot.open("pot.txt");
    vFileOutEfx.open("efx.txt");
    vFileOutNud.open("nud.txt");
    vFileOutEld.open("eld.txt");
    
    G4int imax = 8192;
    G4double vXposition = 0.;
    G4double vXpositionConstant = 1. * vXtalPhysLattice->ComputeInterplanarPeriod();
    
    for(G4int i = 0;i<imax;i++){
        vXposition = double(i) / double(imax) * vXpositionConstant;
        vFileOutPot << vXposition / angstrom << " " << (fPotentialEnergy->ComputeValue(G4ThreeVector(vXposition,0.,0.),vXtalPhysLattice)).x() / eV << std::endl;
        vFileOutEfx << vXposition / angstrom << " " << (fElectricField->ComputeValue(G4ThreeVector(vXposition,0.,0.),vXtalPhysLattice)).x() / eV * angstrom << std::endl;
        vFileOutNud << vXposition / angstrom << " " << (fNucleiDensity->ComputeValue(G4ThreeVector(vXposition,0.,0.),vXtalPhysLattice)).x() * angstrom << std::endl;
        vFileOutEld << vXposition / angstrom << " " << (fElectronDensity->ComputeValue(G4ThreeVector(vXposition,0.,0.),vXtalPhysLattice)).x() * angstrom << std::endl;
    }
    
    vFileOutPot.close();
    vFileOutEfx.close();
    vFileOutNud.close();
    vFileOutEld.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

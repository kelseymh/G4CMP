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

#include "ChannelingParticleUserInfo.hh"

ProcessChanneling::ProcessChanneling(const G4String& aName):G4VDiscreteProcess(aName){
    fLatticeManager = XLatticeManager3::GetXLatticeManager();
    
    G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
    G4cout<<"\n ProcessChanneling::Constructor: Geometry surface tolerance is: " << kCarTolerance / mm << " mm"<<std::endl;
    if(verboseLevel>1) {
        G4cout << GetProcessName() << " is created "<< G4endl;
    }
        
    fFileOut.open("channelingNC.txt");
    fFileOut << "index,posin,angin,depth,pos,ang,dens,tr_en,ndch" << std::endl;
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

XVCrystalIntegratedDensity* ProcessChanneling::GetIntegratedDensityNuclei(){
    return fIntegratedDensityNuclei;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ProcessChanneling::SetIntegratedDensityNuclei(XVCrystalIntegratedDensity* vIntegratedDensity){
    fIntegratedDensityNuclei = vIntegratedDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XVCrystalIntegratedDensity* ProcessChanneling::GetIntegratedDensityElectron(){
    return fIntegratedDensityElectron;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ProcessChanneling::SetIntegratedDensityElectron(XVCrystalIntegratedDensity* vIntegratedDensity){
    fIntegratedDensityElectron = vIntegratedDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ProcessChanneling::UpdatePositionMomentumDensity(const G4Track& aTrack){    
    
    if(!fIntegratedDensityNuclei->HasBeenInitialized(GetXPhysicalLattice(aTrack)) ||
       !fIntegratedDensityElectron->HasBeenInitialized(GetXPhysicalLattice(aTrack)) ){
        ComputeCrystalCharacteristicForChanneling(aTrack);
        G4cout << "ChannelingProcess::UpdatePositionMomentumDensity::fIntegratedDensityNuclei->Initialized" << std::endl;
        G4cout << "ChannelingProcess::UpdatePositionMomentumDensity::fIntegratedDensityElectron->Initialized" << std::endl;
    }

    if(!GetInfo(aTrack)->GetChanneling()){
        
        if(GetInfo(aTrack)->GetPositionChanneledInitial().x() == DBL_MAX){
            G4double vXposition = G4UniformRand() * GetXPhysicalLattice(aTrack)->ComputeInterplanarPeriod();
            GetInfo(aTrack)->SetPositionChanneled(G4ThreeVector(vXposition,0.,0.));
            GetInfo(aTrack)->SetPositionChanneledInitial(GetInfo(aTrack)->GetPositionChanneled());
        }
        else{
            GetInfo(aTrack)->SetPositionChanneled(ComputeNewPosition(aTrack));
        }

        GetInfo(aTrack)->SetMomentumChanneled(fLatticeManager->GetXPhysicalLattice(GetVolume(aTrack))->ProjectVectorFromWorldToLattice(aTrack.GetMomentum()));

        if(GetInfo(aTrack)->GetMomentumChanneledInitial().x() == DBL_MAX){
            GetInfo(aTrack)->SetMomentumChanneledInitial(GetInfo(aTrack)->GetMomentumChanneled());
        }
    }
    else{
        G4ThreeVector vMomentum = GetInfo(aTrack)->GetMomentumChanneled();
        vMomentum += fLatticeManager->GetXPhysicalLattice(GetVolume(aTrack))->ProjectVectorFromWorldToLattice(aTrack.GetMomentum());
        GetInfo(aTrack)->SetMomentumChanneled(vMomentum);
    }
    
    GetInfo(aTrack)->SetNucleiDensityPreviousStep(GetInfo(aTrack)->GetNucleiDensity());
    GetInfo(aTrack)->SetNucleiDensity(fIntegratedDensityNuclei->GetValue(ComputeTransverseEnergy(aTrack).x(),GetXPhysicalLattice(aTrack)));
    
    GetInfo(aTrack)->SetElectronDensityPreviousStep(GetInfo(aTrack)->GetElectronDensity());
    GetInfo(aTrack)->SetElectronDensity(fIntegratedDensityElectron->GetValue(ComputeTransverseEnergy(aTrack).x(),GetXPhysicalLattice(aTrack)));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool ProcessChanneling::IsInChanneling(const G4Track& aTrack){
    //----------------------------------------
    // check if the particle momentum
    // transverse to the (h,k,l) plane
    // is small enough to permit channeling
    //----------------------------------------
    
    
    UpdatePositionMomentumDensity(aTrack);
    if( ComputeTransverseEnergy(aTrack).x() <= ComputeChannelingCriticalEnergy(aTrack) ){
        GetInfo(aTrack)->SetChanneling(true);
        return true;
    }
    return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double ProcessChanneling::ComputeChannelingCriticalEnergy(const G4Track& aTrack){
    //----------------------------------------
    // compute the critical energy
    // for channeling
    //----------------------------------------
        
    G4double vCriticalEnergy = 0.;
    vCriticalEnergy += fPotentialEnergy->GetMaximum(GetXPhysicalLattice(aTrack));

    return vCriticalEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector ProcessChanneling::ComputeTransverseEnergy(const G4Track& aTrack){
    //----------------------------------------
    // compute the particle transverse energy
    // in the crystal reference system
    //----------------------------------------
            
    G4ThreeVector vTransverseEnergy = G4ThreeVector(0.,0.,0.);
    
    G4double vTotalEnergy = aTrack.GetStep()->GetPostStepPoint()->GetTotalEnergy();
    
    vTransverseEnergy += fPotentialEnergy->ComputeValue(GetInfo(aTrack)->GetPositionChanneled(),GetXPhysicalLattice(aTrack));
    
    vTransverseEnergy += G4ThreeVector( 0., 0., (  pow( GetInfo(aTrack)->GetMomentumChanneled().z(), 2. )  / vTotalEnergy ) );
    
    vTransverseEnergy += G4ThreeVector( (  pow( GetInfo(aTrack)->GetMomentumChanneled().x(), 2. ) / vTotalEnergy ), 0., 0. );
    
    return vTransverseEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector ProcessChanneling::ComputeChannelingOutgoingMomentum(const G4Track& aTrack){
    
    G4double vTotalEnergy = aTrack.GetStep()->GetPostStepPoint()->GetTotalEnergy();
    
    G4double vPhi = G4UniformRand() * GetInfo(aTrack)->GetMomentumChanneled().x()/vTotalEnergy;
    
    G4double vTheta = G4UniformRand() * GetInfo(aTrack)->GetMomentumChanneled().z()/vTotalEnergy;
    
    G4ThreeVector vNewMomentum = G4ThreeVector(0.,1.,0.);
    
    vNewMomentum.rotate(G4ThreeVector(0,0,1),vPhi).rotate(G4ThreeVector(1,0,0),vTheta);
    
    vNewMomentum = GetXPhysicalLattice(aTrack)->ProjectVectorFromLatticeToWorld(vNewMomentum);

    vNewMomentum = vNewMomentum.unit();

    return vNewMomentum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double ProcessChanneling::GetChannelingMeanFreePath(const G4Track& aTrack){
    //----------------------------------------
    // return the channeling MFP
    //----------------------------------------
    
    G4double vFactor = 10.;
    
    G4double vMFP = vFactor * ComputeOscillationPeriod(aTrack);
    if(GetInfo(aTrack)->GetNucleiDensity() != 0){
        vMFP /= GetInfo(aTrack)->GetNucleiDensity();
    }
    else{
        vMFP * 100.;
    }

    return vMFP;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double ProcessChanneling::GetMeanFreePath(const G4Track& aTrack,
                                            G4double previousStepSize,
                                            G4ForceCondition* condition){
    
    //----------------------------------------
    // the condition is forced to check if
    // the volume has a lattice at each step.
    // if it hasn't, return DBL_MAX
    //----------------------------------------
    
    *condition = Forced;
        
    if(fLatticeManager->HasLattice(GetVolume(aTrack))){
        return GetChannelingMeanFreePath(aTrack);
    }
    else{
        return DBL_MAX;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* ProcessChanneling::PostStepDoIt(const G4Track& aTrack,
                                                   const G4Step& aStep){
    
    //----------------------------------------
    // check if the volume has a lattice
    // and if the particle is in channeling.
    // If it is so, the particle is forced
    // to follow the channeling plane
    // direction. If the particle has
    // dechanneled or exited the crystal
    // the outgoing angle is evaluated
    //----------------------------------------
    
    aParticleChange.Initialize(aTrack);
            
    bool bIsInChanneling = false;
    
    if(fLatticeManager->HasLattice(GetVolume(aTrack)))
    {
        bIsInChanneling = IsInChanneling(aTrack);
        if(bIsInChanneling){
            aParticleChange.ProposeMomentumDirection(GetXPhysicalLattice(aTrack)->GetLatticeDirection().unit());
        }
    }
    
    if(!bIsInChanneling)
    {
        if(GetInfo(aTrack)->GetChanneling()){            
            aParticleChange.ProposeMomentumDirection(ComputeChannelingOutgoingMomentum(aTrack));
            
            fFileOut << G4int(aTrack.GetTrackID()) << "," << GetInfo(aTrack)->GetPositionChanneledInitial().x()/angstrom << "," << GetInfo(aTrack)->GetMomentumChanneledInitial().x()/aTrack.GetStep()->GetPostStepPoint()->GetTotalEnergy()*1.E6   << "," << aTrack.GetStep()->GetPostStepPoint()->GetPosition().y()/micrometer << "," << GetInfo(aTrack)->GetPositionChanneled().x()/angstrom << "," << GetInfo(aTrack)->GetMomentumChanneled().x()/aTrack.GetStep()->GetPostStepPoint()->GetTotalEnergy()*1.E6 << "," << GetInfo(aTrack)->GetNucleiDensity() << "," <<  ComputeTransverseEnergy(aTrack).x()/eV << "," << GetInfo(aTrack)->GetNumberOfDechanneling()<< std::endl;
            
            GetInfo(aTrack)->IncreaseNumberOfDechanneling();
        }
        GetInfo(aTrack)->SetChanneling(false);
        GetInfo(aTrack)->SetNucleiDensity(1.);
    }
    
    return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4bool ProcessChanneling::IsApplicable(const G4ParticleDefinition& aPD){
    return(aPD.GetPDGCharge()>0.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void ProcessChanneling::BuildPhysicsTable(const G4ParticleDefinition&){
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double ProcessChanneling::ComputeChannelingCriticalAngle(const G4Track& aTrack){
    //----------------------------------------
    // compute the critical angle
    // for chenneling
    //----------------------------------------
    
    G4double vTotalEnergy = aTrack.GetStep()->GetPostStepPoint()->GetTotalEnergy();

    G4double vCriticalAngle = pow( 2.0 * fabs( ComputeChannelingCriticalEnergy(aTrack) / vTotalEnergy ) , 0.5);
    
    return vCriticalAngle;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector ProcessChanneling::ComputeNewPosition(const G4Track& aTrack){
    return GetInfo(aTrack)->GetPositionChanneledInitial();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double ProcessChanneling::ComputeOscillationPeriod(const G4Track& aTrack){
    //----------------------------------------
    // compute the particle oscillation
    // period in the crystal channel
    //----------------------------------------

    G4double vInterplanarPeriod = GetXPhysicalLattice(aTrack)->ComputeInterplanarPeriod();
    G4double vOscillationPeriod = M_PI * vInterplanarPeriod / ComputeChannelingCriticalAngle(aTrack);
    return vOscillationPeriod;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XPhysicalLattice* ProcessChanneling::GetXPhysicalLattice(const G4Track& aTrack){
    return fLatticeManager->GetXPhysicalLattice(aTrack.GetStep()->GetPostStepPoint()->GetPhysicalVolume());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* ProcessChanneling::GetVolume(const G4Track& aTrack){
    return aTrack.GetStep()->GetPostStepPoint()->GetPhysicalVolume();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ChannelingParticleUserInfo* ProcessChanneling::GetInfo(const G4Track& aTrack){
    return (ChannelingParticleUserInfo*) aTrack.GetUserInformation();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ProcessChanneling::ComputeCrystalCharacteristicForChanneling(const G4Track& aTrack){
    
    fIntegratedDensityNuclei->SetXPhysicalLattice(GetXPhysicalLattice(aTrack));
    fIntegratedDensityNuclei->InitializeTable();
    
    fIntegratedDensityElectron->SetXPhysicalLattice(GetXPhysicalLattice(aTrack));
    fIntegratedDensityElectron->InitializeTable();
    
    PrintCrystalCharacteristicsOnFiles(aTrack);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ProcessChanneling::PrintCrystalCharacteristicsOnFiles(const G4Track& aTrack){
    char* filename;
    fIntegratedDensityNuclei->PrintOnFile(filename="dens_nuclei.txt",GetXPhysicalLattice(aTrack));
    fIntegratedDensityElectron->PrintOnFile(filename="dens_electrons.txt",GetXPhysicalLattice(aTrack));
    fPotentialEnergy->PrintOnFile(filename="pot.txt",GetXPhysicalLattice(aTrack),eV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

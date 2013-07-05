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


#include "G4Box.hh"

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
    
    fFileName = "";
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

void ProcessChanneling::SetFileName(const G4String& vFilename){
    
    if(fFileOut && vFilename && vFilename!=fFileName){
        fFileOut.close();
        
        fFileName = vFilename;
        
        G4String filename;
        fFileOut.open(filename = fFileName + "_CH.txt");
        fFileOut << "index,posin,angin,depth,pos,ang,dens,tr_en,ndch" << std::endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4String ProcessChanneling::GetFileName(){
    return fFileName;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ProcessChanneling::PrintChannelingPropertiesOnFile(const G4Track& aTrack){
    fFileOut << G4int(aTrack.GetTrackID()) << "," << GetInfo(aTrack)->GetPositionChanneledInitial().x()/angstrom << "," << GetInfo(aTrack)->GetMomentumChanneledInitial().x()/aTrack.GetStep()->GetPostStepPoint()->GetTotalEnergy()*1.E6   << "," << aTrack.GetStep()->GetPostStepPoint()->GetPosition().y()/micrometer << "," << GetInfo(aTrack)->GetPositionChanneled().x()/angstrom << "," << GetInfo(aTrack)->GetMomentumChanneled().x()/aTrack.GetStep()->GetPostStepPoint()->GetTotalEnergy()*1.E6 << "," << GetInfo(aTrack)->GetNucleiDensity() << "," <<  ComputeTransverseEnergy(aTrack).x()/eV << "," << GetInfo(aTrack)->GetNumberOfDechanneling()<< std::endl;
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

XCrystalIntegratedDensityHub* ProcessChanneling::GetIntegratedDensity(){
    return fIntegratedDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ProcessChanneling::SetIntegratedDensity(XCrystalIntegratedDensityHub* vIntegratedDensity){
    fIntegratedDensity = vIntegratedDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ProcessChanneling::UpdatePositionMomentumDensity(const G4Track& aTrack){
    
    if(!fIntegratedDensity->HasBeenInitialized(GetXPhysicalLattice(aTrack))){
        ComputeCrystalCharacteristic(aTrack);
        G4cout << "ChannelingProcess::UpdatePositionMomentumDensity::fIntegratedDensity->Initialized" << std::endl;
    }
    
    UpdatePosition(aTrack);
    UpdateMomentum(aTrack);
    UpdateDensity(aTrack);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ProcessChanneling::UpdatePosition(const G4Track& aTrack){
    if(!GetInfo(aTrack)->HasBeenUnderCoherentEffect()){
        if(GetInfo(aTrack)->GetPositionChanneledInitial().x() == DBL_MAX){
            G4double vXposition = G4UniformRand() * GetXPhysicalLattice(aTrack)->ComputeInterplanarPeriod();
            GetInfo(aTrack)->SetPositionChanneled(G4ThreeVector(vXposition,0.,0.));
            GetInfo(aTrack)->SetPositionChanneledInitial(GetInfo(aTrack)->GetPositionChanneled());
        }
        else{
            GetInfo(aTrack)->SetPositionChanneled(ComputePostStepPositionInTheChannel(aTrack));
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ProcessChanneling::UpdateMomentum(const G4Track& aTrack){
    if(GetInfo(aTrack)->HasBeenUnderCoherentEffect() == 1){
        G4ThreeVector vMomentum = GetInfo(aTrack)->GetMomentumChanneled();
        
        G4ThreeVector vPosition = ComputePositionInTheCrystal(aTrack.GetStep()->GetPreStepPoint(),aTrack);
        
        vMomentum += GetXPhysicalLattice(aTrack)->ProjectMomentumVectorFromWorldToLattice(aTrack.GetMomentum(),vPosition);
        
        GetInfo(aTrack)->SetMomentumChanneled(vMomentum);
    }
    else{
        
        G4ThreeVector vPosition = ComputePositionInTheCrystal(aTrack.GetStep()->GetPreStepPoint(),aTrack);
        // we take the PREVIOUS step point to compare, otherwise the momentum is not computed correctly
        
        G4ThreeVector vMomentumChanneling = GetXPhysicalLattice(aTrack)->ProjectMomentumVectorFromWorldToLattice(aTrack.GetMomentum(),vPosition);
        
        GetInfo(aTrack)->SetMomentumChanneled(vMomentumChanneling);
        
        if(GetInfo(aTrack)->GetMomentumChanneledInitial().x() == DBL_MAX){
            GetInfo(aTrack)->SetMomentumChanneledInitial(GetInfo(aTrack)->GetMomentumChanneled());
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ProcessChanneling::UpdateDensity(const G4Track& aTrack){
    GetInfo(aTrack)->SetNucleiDensityPreviousStep(GetInfo(aTrack)->GetNucleiDensity());
    GetInfo(aTrack)->SetElectronDensityPreviousStep(GetInfo(aTrack)->GetElectronDensity());
    
    if(GetXPhysicalLattice(aTrack)->IsBent() == false){
        GetInfo(aTrack)->SetNucleiDensity(fIntegratedDensity->GetIntegratedDensityNuclei(ComputeTransverseEnergy(aTrack).x(),GetXPhysicalLattice(aTrack),GetParticleDefinition(aTrack)->GetPDGCharge()));
        GetInfo(aTrack)->SetElectronDensity(fIntegratedDensity->GetIntegratedDensityElectron(ComputeTransverseEnergy(aTrack).x(),GetXPhysicalLattice(aTrack),GetParticleDefinition(aTrack)->GetPDGCharge()));
    }
    else{
        GetInfo(aTrack)->SetNucleiDensity(fIntegratedDensity->GetIntegratedDensityNuclei(ComputeTransverseEnergyBent(aTrack).x(),GetXPhysicalLattice(aTrack),GetParticleDefinition(aTrack)->GetPDGCharge()));
        GetInfo(aTrack)->SetElectronDensity(fIntegratedDensity->GetIntegratedDensityElectron(ComputeTransverseEnergyBent(aTrack).x(),GetXPhysicalLattice(aTrack),GetParticleDefinition(aTrack)->GetPDGCharge()));
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool ProcessChanneling::IsUnderCoherentEffect(const G4Track& aTrack){
    //----------------------------------------
    // check if the particle momentum
    // transverse to the (h,k,l) plane
    // is small enough to permit channeling
    //----------------------------------------
    
    
    UpdatePositionMomentumDensity(aTrack);
    
    G4double vCriticalEnergy = 0.;
    G4double vTransverseEnergy = 0.;
    
    if(GetXPhysicalLattice(aTrack)->IsBent() == false){
        vCriticalEnergy = ComputeCriticalEnergy(aTrack);
        vTransverseEnergy = ComputeTransverseEnergy(aTrack).x();
        if( vTransverseEnergy <= vCriticalEnergy ){
            GetInfo(aTrack)->SetCoherentEffect(1);
            return true;
        }
    }
    else{
        vCriticalEnergy = ComputeCriticalEnergyBent(aTrack);
        vTransverseEnergy = ComputeTransverseEnergyBent(aTrack).x();
        if( vTransverseEnergy <=  vCriticalEnergy){
            vCriticalEnergy = ComputeCriticalEnergy(aTrack);
           if( vTransverseEnergy <= vCriticalEnergy ){
                GetInfo(aTrack)->SetCoherentEffect(1);
            }
            else{
                GetInfo(aTrack)->SetCoherentEffect(2);
            }
            return true;
        }
    }
    
    return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double ProcessChanneling::ComputeCriticalEnergy(const G4Track& aTrack){
    //----------------------------------------
    // compute the critical energy
    // for channeling
    //----------------------------------------
    
    G4double vCriticalEnergy = 0.;
    vCriticalEnergy += fPotentialEnergy->GetMaximum(GetXPhysicalLattice(aTrack));
    
    return vCriticalEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double ProcessChanneling::ComputeCriticalEnergyBent(const G4Track& aTrack){
    //----------------------------------------
    // compute the critical energy
    // for channeling
    //----------------------------------------
    
    
    G4double vTotalEnergy = aTrack.GetStep()->GetPostStepPoint()->GetTotalEnergy();
    G4double vInterplanarPeriod = GetXPhysicalLattice(aTrack)->ComputeInterplanarPeriod();
    
    G4double vCriticalEnergy = ComputeCriticalEnergy(aTrack);
    vCriticalEnergy += vTotalEnergy * vInterplanarPeriod / fabs(GetXPhysicalLattice(aTrack)->GetCurvatureRadius().x());
    
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
    
    G4double vParticleCharge = GetParticleDefinition(aTrack)->GetPDGCharge();
    
    vTransverseEnergy += vParticleCharge * fPotentialEnergy->ComputeEC(GetInfo(aTrack)->GetPositionChanneled(),GetXPhysicalLattice(aTrack));
    
    vTransverseEnergy += G4ThreeVector( 0., 0., (  pow( GetInfo(aTrack)->GetMomentumChanneled().z(), 2. )  / vTotalEnergy ) );
    
    vTransverseEnergy += G4ThreeVector( ( pow( GetInfo(aTrack)->GetMomentumChanneled().x(), 2. ) / vTotalEnergy ), 0., 0. );
    
    return vTransverseEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector ProcessChanneling::ComputeTransverseEnergyBent(const G4Track& aTrack){
    //----------------------------------------
    // compute the particle transverse energy
    // in the crystal reference system
    //----------------------------------------
    
    G4ThreeVector vTransverseEnergy = ComputeTransverseEnergy(aTrack);
    
    G4double vTotalEnergy = aTrack.GetStep()->GetPostStepPoint()->GetTotalEnergy();
    
    G4double vTransverseEnergyX = vTransverseEnergy.x() + vTotalEnergy * GetInfo(aTrack)->GetPositionChanneled().x() / fabs(GetXPhysicalLattice(aTrack)->GetCurvatureRadius().x());
    
    vTransverseEnergy.setX(vTransverseEnergyX);
    
    return vTransverseEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector ProcessChanneling::ComputeChannelingOutgoingMomentum(const G4Track& aTrack){
    
    G4double vTotalEnergy = aTrack.GetStep()->GetPostStepPoint()->GetTotalEnergy();
    
    G4double vPhi = 2. * ( G4UniformRand() - 0.5) * GetInfo(aTrack)->GetMomentumChanneled().x() / vTotalEnergy;
    
    G4double vTheta = 2. * ( G4UniformRand() - 0.5) * GetInfo(aTrack)->GetMomentumChanneled().z() / vTotalEnergy;
    
    G4ThreeVector vNewMomentum = G4ThreeVector(0.,1.,0.);
    
    vNewMomentum.rotate(G4ThreeVector(0,0,1),vPhi).rotate(G4ThreeVector(1,0,0),vTheta);
    
    G4ThreeVector vPosition = ComputePositionInTheCrystal(aTrack.GetStep()->GetPostStepPoint(),aTrack);
    
    vNewMomentum = GetXPhysicalLattice(aTrack)->ProjectMomentumVectorFromLatticeToWorld(vNewMomentum,vPosition);
    
    vNewMomentum = vNewMomentum.unit();
    
    return vNewMomentum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector ProcessChanneling::ComputeVolumeReflectionOutgoingMomentum(const G4Track& aTrack){
    
    G4double vTotalEnergy = aTrack.GetStep()->GetPostStepPoint()->GetTotalEnergy();
    
    G4double vPhi = - 2. * GetInfo(aTrack)->GetMomentumChanneled().x() / vTotalEnergy;
    
    G4ThreeVector vNewMomentum = G4ThreeVector(0.,1.,0.);
    
    vNewMomentum.rotate(G4ThreeVector(0,0,1),vPhi);
    
    G4ThreeVector vPosition = ComputePositionInTheCrystal(aTrack.GetStep()->GetPostStepPoint(),aTrack);
    
    vNewMomentum = GetXPhysicalLattice(aTrack)->ProjectMomentumVectorFromLatticeToWorld(vNewMomentum,vPosition);
    
    vNewMomentum = vNewMomentum.unit();
    
    return vNewMomentum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector ProcessChanneling::ComputePositionInTheCrystal(G4StepPoint* vStep,const G4Track& aTrack){
    G4VPhysicalVolume* vVolumePre = aTrack.GetStep()->GetPreStepPoint()->GetPhysicalVolume();
    G4VPhysicalVolume* vVolumePost = aTrack.GetStep()->GetPostStepPoint()->GetPhysicalVolume();
    
    G4ThreeVector vLocalPos;
    
    if(!fLatticeManager->HasLattice(vVolumePre)){
        vLocalPos = G4ThreeVector(0.,0.,0.);
    }
    
    if(!fLatticeManager->HasLattice(vVolumePost)){
        G4Box* vXtalSolid = (G4Box*) aTrack.GetStep()->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetSolid();
        vLocalPos = G4ThreeVector(0.,vXtalSolid->GetYHalfLength()*2.,0.);
        return vLocalPos;
    }
    
    G4TouchableHistory* theTouchable = (G4TouchableHistory*)(vStep->GetTouchable());
    G4ThreeVector vWorldPos = vStep->GetPosition();
    vLocalPos = theTouchable->GetHistory()->GetTopTransform().TransformPoint(vWorldPos);
    G4Box* vXtalSolid = (G4Box*) vStep->GetPhysicalVolume()->GetLogicalVolume()->GetSolid();
    vLocalPos += G4ThreeVector(vXtalSolid->GetXHalfLength(),vXtalSolid->GetYHalfLength(),vXtalSolid->GetZHalfLength());
    
    return vLocalPos;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double ProcessChanneling::GetChannelingMeanFreePath(const G4Track& aTrack){
    //----------------------------------------
    // return the channeling MFP
    //----------------------------------------
    
    G4double vFactor = 10.;
    
    G4double vMFP = vFactor * ComputeOscillationPeriod(aTrack);
    
    if(GetInfo(aTrack)->GetNucleiDensity() > 0.1){
        vMFP /= GetInfo(aTrack)->GetNucleiDensity();
    }
    else{
        vMFP * 10.;
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
    // dechanneled or exited the crystal,
    // the outgoing angle is evaluated
    //----------------------------------------
    
    aParticleChange.Initialize(aTrack);
    
    bool bIsUnderCoherentEffect = false;
    
    if(GetInfo(aTrack)->HasBeenUnderCoherentEffect() == 0 ||
       !fLatticeManager->HasLattice(GetVolume(aTrack))){
        GetInfo(aTrack)->SetNucleiDensityPreviousStep(1.);
        GetInfo(aTrack)->SetElectronDensityPreviousStep(1.);
        
        GetInfo(aTrack)->SetNucleiDensity(1.);
        GetInfo(aTrack)->SetElectronDensity(1.);
    }
    
    if(fLatticeManager->HasLattice(GetVolume(aTrack))){
        bIsUnderCoherentEffect = IsUnderCoherentEffect(aTrack);
        if(bIsUnderCoherentEffect){
            if(GetInfo(aTrack)->HasBeenUnderCoherentEffect() == 1){
                G4ThreeVector vPosition = ComputePositionInTheCrystal(aTrack.GetStep()->GetPostStepPoint(),aTrack);
                aParticleChange.ProposeMomentumDirection(GetXPhysicalLattice(aTrack)->GetLatticeDirection(vPosition).unit());
            }
            else if(GetInfo(aTrack)->HasBeenUnderCoherentEffect() == 2){
                aParticleChange.ProposeMomentumDirection(ComputeVolumeReflectionOutgoingMomentum(aTrack));
            }
        }
    }
    
    if(!bIsUnderCoherentEffect){
        if(GetInfo(aTrack)->HasBeenUnderCoherentEffect() == 1){
            aParticleChange.ProposeMomentumDirection(ComputeChannelingOutgoingMomentum(aTrack));
            
            if(fFileOut){
                PrintChannelingPropertiesOnFile(aTrack);
            }
            GetInfo(aTrack)->IncreaseNumberOfDechanneling();
        }
        if(GetInfo(aTrack)->HasBeenUnderCoherentEffect() == 2){
            if(fFileOut){
                PrintChannelingPropertiesOnFile(aTrack);
            }
        }

        
        GetInfo(aTrack)->SetCoherentEffect(false);
        
        GetInfo(aTrack)->SetNucleiDensity(1.);
        GetInfo(aTrack)->SetElectronDensity(1.);
    }
    
    return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4bool ProcessChanneling::IsApplicable(const G4ParticleDefinition& aPD){
    return(aPD.GetPDGCharge() != 0.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void ProcessChanneling::BuildPhysicsTable(const G4ParticleDefinition&){
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double ProcessChanneling::ComputeCriticalAngle(const G4Track& aTrack){
    //----------------------------------------
    // compute the critical angle
    // for chenneling
    //----------------------------------------
    
    G4double vTotalEnergy = aTrack.GetStep()->GetPostStepPoint()->GetTotalEnergy();
    
    G4double vCriticalAngle = pow( 2.0 * fabs( ComputeCriticalEnergy(aTrack) / vTotalEnergy ) , 0.5);
    
    return vCriticalAngle;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector ProcessChanneling::ComputePostStepPositionInTheChannel(const G4Track& aTrack){
    return GetInfo(aTrack)->GetPositionChanneledInitial();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double ProcessChanneling::ComputeOscillationPeriod(const G4Track& aTrack){
    //----------------------------------------
    // compute the particle oscillation
    // period in the crystal channel
    //----------------------------------------
    
    G4double vInterplanarPeriod = GetXPhysicalLattice(aTrack)->ComputeInterplanarPeriod();
    G4double vOscillationPeriod = M_PI * vInterplanarPeriod / ComputeCriticalAngle(aTrack);
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

G4ParticleDefinition* ProcessChanneling::GetParticleDefinition(const G4Track& aTrack){
    return const_cast<G4ParticleDefinition*>(aTrack.GetParticleDefinition());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ProcessChanneling::ComputeCrystalCharacteristic(const G4Track& aTrack){
    
    fIntegratedDensity->SetXPhysicalLattice(GetXPhysicalLattice(aTrack));
    fIntegratedDensity->InitializeTables();
    
    if(fFileName != "") {
        PrintCrystalCharacteristicsOnFiles(aTrack);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ProcessChanneling::PrintCrystalCharacteristicsOnFiles(const G4Track& aTrack){
    
    G4String filename;

    fIntegratedDensity->PrintOnFiles(filename=fFileName);
    fPotentialEnergy->PrintOnFile(filename=fFileName + "_pot.txt",GetXPhysicalLattice(aTrack),eV);
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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
        fFileOut << "index,posin,angin,depth,pos,ang,angout,dens,tr_en,ndch,cheff" << std::endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4String ProcessChanneling::GetFileName(){
    return fFileName;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ProcessChanneling::PrintChannelingPropertiesOnFile(const G4Track& aTrack){
    if(fFileOut){
        fFileOut << G4int(aTrack.GetTrackID()) << ","
        << GetInfo(aTrack)->GetPositionChanneledInitial().x() / angstrom << ","
        << GetInfo(aTrack)->GetMomentumChanneledInitial().x() / MeV << ","
        << aTrack.GetStep()->GetPostStepPoint()->GetPosition().y() / micrometer << ","
        << GetInfo(aTrack)->GetPositionChanneled().x() / angstrom << ","
        << GetInfo(aTrack)->GetMomentumChanneled().x() / MeV << ","
        << aTrack.GetMomentum().x() / MeV << ","
        << GetInfo(aTrack)->GetNucleiDensity() << ","
        << ComputeTransverseEnergy(aTrack).x() / eV << ","
        << GetInfo(aTrack)->GetNumberOfDechanneling() << ","
        << GetInfo(aTrack)->HasBeenUnderCoherentEffect() << std::endl;
    }
    GetInfo(aTrack)->IncreaseNumberOfDechanneling();
    
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

XVCrystalCharacteristic* ProcessChanneling::GetElectricField(){
    return fElectricField;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ProcessChanneling::SetElectricField(XVCrystalCharacteristic* vElectricField){
    fElectricField = vElectricField;
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
    
    if(fIntegratedDensity->HasBeenInitialized(GetXPhysicalLattice(aTrack)) == false){
        ComputeCrystalCharacteristic(aTrack);
        G4cout << "ChannelingProcess::UpdatePositionMomentumDensity::fIntegratedDensity->Initialized" << std::endl;
    }
    
    UpdatePosition(aTrack);
    UpdateMomentum(aTrack);
    UpdateDensity(aTrack);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ProcessChanneling::UpdatePosition(const G4Track& aTrack){
    if(GetInfo(aTrack)->GetPositionChanneledInitial().x() == DBL_MAX){
        G4double vXposition = G4UniformRand() * GetXPhysicalLattice(aTrack)->ComputeInterplanarPeriod();
        
        GetInfo(aTrack)->SetPositionChanneled(G4ThreeVector(0.,0.,0.));
        
        GetInfo(aTrack)->SetPositionChanneledInitial(G4ThreeVector(vXposition,0.,0.));
    }
    else{
        if(GetInfo(aTrack)->HasBeenUnderCoherentEffect() == false){
            G4ThreeVector vDifferencePosition = ComputePositionInTheCrystal(aTrack.GetStep()->GetPostStepPoint(),aTrack) - ComputePositionInTheCrystal(aTrack.GetStep()->GetPreStepPoint(),aTrack);
            
            GetInfo(aTrack)->SetPositionChanneled(vDifferencePosition);
        }
        else{
            GetInfo(aTrack)->SetPositionChanneled(G4ThreeVector(0.,0.,0.));
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ProcessChanneling::UpdateMomentum(const G4Track& aTrack){
    if(GetInfo(aTrack)->GetMomentumChanneledInitial().x() == DBL_MAX){
        G4ThreeVector vPosition = ComputePositionInTheCrystal(aTrack.GetStep()->GetPreStepPoint(),aTrack);
        // we take the PREVIOUS step point to compare, otherwise the momentum is not computed correctly
        
        G4ThreeVector vMomentum = GetXPhysicalLattice(aTrack)->ProjectMomentumVectorFromWorldToLattice(aTrack.GetMomentum(),vPosition);
        
        GetInfo(aTrack)->SetMomentumChanneled(vMomentum);
        
        GetInfo(aTrack)->SetEnergyChanneled(0.);

        GetInfo(aTrack)->SetEnergyChanneledInitial(0.);

        if(GetInfo(aTrack)->GetMomentumChanneledInitial().x() == DBL_MAX){
            GetInfo(aTrack)->SetMomentumChanneledInitial(GetInfo(aTrack)->GetMomentumChanneled());
        }
    }
    else{
        G4ThreeVector vMomentum = GetInfo(aTrack)->GetMomentumChanneled();
        
        G4ThreeVector vPosition = ComputePositionInTheCrystal(aTrack.GetStep()->GetPreStepPoint(),aTrack);
        
        vMomentum += GetXPhysicalLattice(aTrack)->ProjectMomentumVectorFromWorldToLattice(aTrack.GetMomentum(),vPosition);
        
        GetInfo(aTrack)->SetMomentumChanneled(vMomentum);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ProcessChanneling::UpdateDensity(const G4Track& aTrack){
    
    G4double vCharge = GetParticleDefinition(aTrack)->GetPDGCharge();
    
    G4double vTransverseEnergy = 0.;
    
    if(GetXPhysicalLattice(aTrack)->IsBent() == false){
        G4double vTransverseEnergy = ComputeTransverseEnergy(aTrack).x();
    }
    else{
        G4double vTransverseEnergy = ComputeTransverseEnergyBent(aTrack).x();
    }
    
    GetInfo(aTrack)->SetNucleiDensity(fIntegratedDensity->GetIntegratedDensityNuclei(vTransverseEnergy,GetXPhysicalLattice(aTrack),vCharge));
    
    GetInfo(aTrack)->SetElectronDensity(fIntegratedDensity->GetIntegratedDensityElectron(vTransverseEnergy,GetXPhysicalLattice(aTrack),vCharge));
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ProcessChanneling::ResetDensity(const G4Track& aTrack){
    GetInfo(aTrack)->SetNucleiDensity(1.);
    GetInfo(aTrack)->SetElectronDensity(1.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool ProcessChanneling::IsUnderCoherentEffect(const G4Track& aTrack){
    //----------------------------------------
    // check if the particle momentum
    // transverse to the (h,k,l) plane
    // is small enough to permit channeling
    //----------------------------------------
    
    
    UpdatePositionMomentumDensity(aTrack);
    
    G4double vCrEnMax = 0.;
    G4double vCrEnMin = 0.;
    G4double vCrEnRange = 0.;
    
    G4double vTransverseEnergy = 0.;
    
    if(GetXPhysicalLattice(aTrack)->IsBent() == false){
        vCrEnMax = ComputeCriticalEnergy(aTrack);
        vCrEnMin = ComputeCriticalEnergyMinimum(aTrack);
        vCrEnRange = vCrEnMax - vCrEnMin;
        vTransverseEnergy = ComputeTransverseEnergy(aTrack).x();
        G4double vPotentialEnergy = ComputePotentialEnergy(aTrack);

        if(vTransverseEnergy <= vCrEnMax &&
           vTransverseEnergy >= vCrEnMin ){
            GetInfo(aTrack)->SetCoherentEffect(1);
            GetInfo(aTrack)->SetEnergyChanneledInitial(GetInfo(aTrack)->GetEnergyChanneled());
            GetInfo(aTrack)->SetEnergyChanneled(vPotentialEnergy);
            return true;
        }
        else{
            return false;
        }
    }
    else{
        vCrEnMax = ComputeCriticalEnergy(aTrack);
        vCrEnMin = ComputeCriticalEnergyMinimumBent(aTrack);
        vCrEnRange = vCrEnMax - vCrEnMin;
        vTransverseEnergy = ComputeTransverseEnergyBent(aTrack).x();
        
        G4double vCrEnMaxBent = ComputeCriticalEnergyBent(aTrack);

        G4double vPotentialEnergy = ComputePotentialEnergyBent(aTrack);

        G4double vRadiusX = GetXPhysicalLattice(aTrack)->GetCurvatureRadius().x();
        
        if(vPotentialEnergy<0)        G4cout << vTransverseEnergy/eV << " " << vPotentialEnergy/eV << " VR1" << G4endl;

        if(vTransverseEnergy <= vCrEnMax &&
           vTransverseEnergy >= vCrEnMin ){
            GetInfo(aTrack)->SetCoherentEffect(1);
            GetInfo(aTrack)->SetEnergyChanneledInitial(GetInfo(aTrack)->GetEnergyChanneled());
            GetInfo(aTrack)->SetEnergyChanneled(vPotentialEnergy);
            return true;
        }
        else if(vPotentialEnergy > GetInfo(aTrack)->GetEnergyChanneled() &&
                GetInfo(aTrack)->GetEnergyChanneledInitial() > GetInfo(aTrack)->GetEnergyChanneled()){
//            G4cout << vTransverseEnergy/eV << " " << vPotentialEnergy/eV << " VR1" << G4endl;
//            G4cout << GetInfo(aTrack)->GetEnergyChanneled()/eV << " " << GetInfo(aTrack)->GetEnergyChanneledInitial()/eV << " VR2" << G4endl;
            GetInfo(aTrack)->SetCoherentEffect(2);
            GetInfo(aTrack)->SetEnergyChanneledInitial(GetInfo(aTrack)->GetEnergyChanneled());
            GetInfo(aTrack)->SetEnergyChanneled(vPotentialEnergy);
            return true;
        }
        else{
            GetInfo(aTrack)->SetCoherentEffect(0);
            GetInfo(aTrack)->SetEnergyChanneledInitial(GetInfo(aTrack)->GetEnergyChanneled());
            GetInfo(aTrack)->SetEnergyChanneled(vPotentialEnergy);
            return false;
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
    
    
    G4double vTotalEnergy = aTrack.GetStep()->GetPreStepPoint()->GetTotalEnergy();
    
    G4double vInterplanarPeriod = GetXPhysicalLattice(aTrack)->ComputeInterplanarPeriod();
    
    G4double vCriticalEnergy = ComputeCriticalEnergy(aTrack);
    
    vCriticalEnergy += vTotalEnergy * vInterplanarPeriod / fabs(GetXPhysicalLattice(aTrack)->GetCurvatureRadius().x());
    
    return vCriticalEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double ProcessChanneling::ComputeCriticalEnergyMinimum(const G4Track& aTrack){
    //----------------------------------------
    // compute the critical energy minimum
    // for channeling
    //----------------------------------------
    
    G4double vCriticalEnergy = 0.;
    
    vCriticalEnergy += fPotentialEnergy->GetMinimum(GetXPhysicalLattice(aTrack));
    
    return vCriticalEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double ProcessChanneling::ComputeCriticalEnergyMinimumBent(const G4Track& aTrack){
    //----------------------------------------
    // compute the critical energy
    // for channeling
    //----------------------------------------
    
    G4double vTotalEnergy = aTrack.GetStep()->GetPreStepPoint()->GetTotalEnergy();
    
    G4double vInterplanarPeriod = GetXPhysicalLattice(aTrack)->ComputeInterplanarPeriod();
    
    
    G4double vCentre = vInterplanarPeriod * 0.5;
    
    vCentre *= (1. - ComputeCriticalRadius(aTrack) / GetXPhysicalLattice(aTrack)->GetCurvatureRadius().x() );
    
    
    G4double vCriticalEnergy = ComputeCriticalEnergyMinimum(aTrack);
    
    vCriticalEnergy += vTotalEnergy * vCentre / GetXPhysicalLattice(aTrack)->GetCurvatureRadius().x();
    
    return vCriticalEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double ProcessChanneling::ComputeCriticalRadius(const G4Track& aTrack){
    //----------------------------------------
    // compute the critical energy
    // for channeling
    //----------------------------------------
    
    G4double vTotalEnergy = aTrack.GetStep()->GetPreStepPoint()->GetTotalEnergy();
    
    G4double vCriticalRadius = vTotalEnergy / fElectricField->GetMaximum(GetXPhysicalLattice(aTrack)) ;
    
    return vCriticalRadius;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector ProcessChanneling::ComputeTransverseEnergy(const G4Track& aTrack){
    //----------------------------------------
    // compute the particle transverse energy
    // in the crystal reference system
    //----------------------------------------
    
    G4ThreeVector vTransverseEnergy = G4ThreeVector(0.,0.,0.);
    
    G4double vTotalEnergy = aTrack.GetStep()->GetPostStepPoint()->GetTotalEnergy();
    
    G4double vCharge = GetParticleDefinition(aTrack)->GetPDGCharge();
    
    G4ThreeVector vMom = GetInfo(aTrack)->GetMomentumChanneled();
    
    vTransverseEnergy += vCharge * fPotentialEnergy->ComputeEC(GetInfo(aTrack)->GetPositionChanneledInitial(),GetXPhysicalLattice(aTrack));
    
    vTransverseEnergy += G4ThreeVector( 0., 0., (  pow( vMom.z(), 2. )  / vTotalEnergy ) );
    
    vTransverseEnergy += G4ThreeVector( ( pow( vMom.x(), 2. ) / vTotalEnergy ), 0., 0. );
    
    return vTransverseEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector ProcessChanneling::ComputeTransverseEnergyBent(const G4Track& aTrack){
    //----------------------------------------
    // compute the particle transverse energy
    // in the crystal reference system
    //----------------------------------------
    
    G4ThreeVector vTransverseEnergy = ComputeTransverseEnergy(aTrack);
    
    G4double vTotalEnergy = aTrack.GetStep()->GetPreStepPoint()->GetTotalEnergy();
    
    G4ThreeVector vPosition = GetInfo(aTrack)->GetPositionChanneled() + GetInfo(aTrack)->GetPositionChanneledInitial();
    
    G4double vTransverseEnergyX = vTransverseEnergy.x() + vTotalEnergy * vPosition.x() / GetXPhysicalLattice(aTrack)->GetCurvatureRadius().x();
    
    vTransverseEnergy.setX(vTransverseEnergyX);
    
    return vTransverseEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double ProcessChanneling::ComputePotentialEnergy(const G4Track& aTrack){
    //----------------------------------------
    // compute the particle transverse energy
    // in the crystal reference system
    //----------------------------------------
    
    G4double vCharge = GetParticleDefinition(aTrack)->GetPDGCharge();
    
    G4double vPotentialEnergy = vCharge * fPotentialEnergy->ComputeEC(GetInfo(aTrack)->GetPositionChanneledInitial(),GetXPhysicalLattice(aTrack)).x();
        
    return vPotentialEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double ProcessChanneling::ComputePotentialEnergyBent(const G4Track& aTrack){
    //----------------------------------------
    // compute the particle transverse energy
    // in the crystal reference system
    //----------------------------------------
    
    G4double vCharge = GetParticleDefinition(aTrack)->GetPDGCharge();

    G4double vPotentialEnergy = vCharge * fPotentialEnergy->ComputeEC(GetInfo(aTrack)->GetPositionChanneledInitial(),GetXPhysicalLattice(aTrack)).x();
    
    G4double vTotalEnergy = aTrack.GetStep()->GetPreStepPoint()->GetTotalEnergy();
    
    G4ThreeVector vPosition = GetInfo(aTrack)->GetPositionChanneled() + GetInfo(aTrack)->GetPositionChanneledInitial();
    
    vPotentialEnergy += vTotalEnergy * vPosition.x() / GetXPhysicalLattice(aTrack)->GetCurvatureRadius().x();
        
    return vPotentialEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector ProcessChanneling::ComputeChannelingOutgoingMomentum(const G4Track& aTrack){
    
    G4StepPoint* vStepPre = aTrack.GetStep()->GetPreStepPoint();
    
    G4StepPoint* vStepPost = aTrack.GetStep()->GetPostStepPoint();
    
    G4double vTotalEnergy = vStepPre->GetTotalEnergy();
    
    G4double vTransverseEnergy = 0.;
    

    vTransverseEnergy = ComputeTransverseEnergy(aTrack).x();
    
    G4double vPhi = 2. * ( G4UniformRand() - 0.5) * fabs(vTransverseEnergy / vTotalEnergy);
    
    G4double vTheta = 2. * ( G4UniformRand() - 0.5) * fabs(vTransverseEnergy / vTotalEnergy);
    
    G4ThreeVector vNewMomentum = G4ThreeVector(0.,1.,0.);
    
    vNewMomentum.rotate(G4ThreeVector(0,0,1),vPhi).rotate(G4ThreeVector(1,0,0),vTheta);
    
    G4ThreeVector vPosition = ComputePositionInTheCrystal(vStepPost,aTrack);
    
    vNewMomentum = GetXPhysicalLattice(aTrack)->ProjectMomentumVectorFromLatticeToWorld(vNewMomentum,vPosition);
    
    vNewMomentum = vNewMomentum.unit();
    
    return vNewMomentum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector ProcessChanneling::ComputeVolumeReflectionOutgoingMomentum(const G4Track& aTrack){
    
    G4StepPoint* vStep = aTrack.GetStep()->GetPreStepPoint();
    
    G4double vPhi = 0.;
    
    if(GetXPhysicalLattice(aTrack)->IsBent()) {
        G4double vRadiusX = GetXPhysicalLattice(aTrack)->GetCurvatureRadius().x();
        
        G4double vTotalEnergy = vStep->GetTotalEnergy();
                
        G4double vTransverseEnergy = fPotentialEnergy->ComputeEC(GetInfo(aTrack)->GetPositionChanneledInitial(),GetXPhysicalLattice(aTrack)).x();

        vTransverseEnergy += vTotalEnergy * GetInfo(aTrack)->GetPositionChanneledInitial().x() / GetXPhysicalLattice(aTrack)->GetCurvatureRadius().x();

        vPhi = + fabs(vRadiusX)/vRadiusX * pow(+ 2. * fabs(vTransverseEnergy) / vTotalEnergy , 0.5);
        
        G4cout << vPhi / CLHEP::radian * 1.E6 << G4endl;
    }
    
    G4ThreeVector vNewMomentum = G4ThreeVector(0.,1.,0.);
    
    vNewMomentum.rotate(G4ThreeVector(0,0,1),vPhi);
    
    G4ThreeVector vPosition = ComputePositionInTheCrystal(vStep,aTrack);
    
    vNewMomentum = GetXPhysicalLattice(aTrack)->ProjectMomentumVectorFromLatticeToWorld(vNewMomentum,vPosition);
    
    vNewMomentum = vNewMomentum.unit();
    
    return vNewMomentum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector ProcessChanneling::ComputePositionInTheCrystal(G4StepPoint* vStep,const G4Track& aTrack){
    
    G4StepPoint* vStepVol = CheckStepPointLatticeForVolume(vStep,aTrack);
    
    G4StepPoint* vStepPos = CheckStepPointLatticeForPosition(vStep,aTrack);
    
    
    
    G4TouchableHistory* theTouchable = (G4TouchableHistory*)(vStepVol->GetTouchable());
    
    G4ThreeVector vWorldPos = vStepPos->GetPosition();
    
    G4ThreeVector vLocalPos = theTouchable->GetHistory()->GetTopTransform().TransformPoint(vWorldPos);
    
    G4Box* vXtalSolid = (G4Box*) vStepVol->GetPhysicalVolume()->GetLogicalVolume()->GetSolid();
    
    vLocalPos += G4ThreeVector(vXtalSolid->GetXHalfLength(),vXtalSolid->GetYHalfLength(),vXtalSolid->GetZHalfLength());
    
    return vLocalPos;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4StepPoint* ProcessChanneling::CheckStepPointLatticeForVolume(G4StepPoint* vStep, const G4Track& aTrack){
    
    G4StepPoint* vStepPre = aTrack.GetStep()->GetPreStepPoint();
    G4StepPoint* vStepPost = aTrack.GetStep()->GetPostStepPoint();
    
    if( fLatticeManager->HasLattice(vStep->GetPhysicalVolume()) ) {
        return vStep;
    }
    else if(fLatticeManager->HasLattice(vStepPost->GetPhysicalVolume()) == false &&
            fLatticeManager->HasLattice(vStepPre->GetPhysicalVolume()) == true &&
            vStep == vStepPost &&
            vStepPost->GetStepStatus() == fGeomBoundary) {
        return vStepPre;
    }
    else if(fLatticeManager->HasLattice(vStepPre->GetPhysicalVolume()) == false &&
            fLatticeManager->HasLattice(vStepPost->GetPhysicalVolume()) == true &&
            vStep == vStepPre &&
            vStepPre->GetStepStatus() == fGeomBoundary) {
        return vStepPost;
    }
    else{
        return vStep;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4StepPoint* ProcessChanneling::CheckStepPointLatticeForPosition(G4StepPoint* vStep, const G4Track& aTrack){
    
    G4StepPoint* vStepPre = aTrack.GetStep()->GetPreStepPoint();
    G4StepPoint* vStepPost = aTrack.GetStep()->GetPostStepPoint();
    
    if( fLatticeManager->HasLattice(vStep->GetPhysicalVolume()) ) {
        return vStep;
    }
    else if(fLatticeManager->HasLattice(vStepPost->GetPhysicalVolume()) == false &&
            fLatticeManager->HasLattice(vStepPre->GetPhysicalVolume()) == true &&
            vStep == vStepPost &&
            vStepPost->GetStepStatus() == fGeomBoundary) {
        return vStepPost;
    }
    else if(fLatticeManager->HasLattice(vStepPre->GetPhysicalVolume()) == false &&
            fLatticeManager->HasLattice(vStepPost->GetPhysicalVolume()) == true &&
            vStep == vStepPre &&
            vStepPre->GetStepStatus() == fGeomBoundary) {
        return vStepPost;
    }
    else{
        return vStep;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double ProcessChanneling::GetChannelingMeanFreePath(const G4Track& aTrack){
    //----------------------------------------
    // return the channeling MFP
    //----------------------------------------
    
    G4double vMFP = 0.;
    
    if(GetInfo(aTrack)->HasBeenUnderCoherentEffect() == 1){
        vMFP = ComputeOscillationPeriod(aTrack) / GetInfo(aTrack)->GetNucleiDensity();
    }
    else{
        vMFP = ComputeOscillationPeriod(aTrack);
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
    
    if(HasLattice(aTrack)){
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
    
    GetInfo(aTrack)->StoreDensityPreviousStep();
    
    G4bool bIsUnderCoherentEffect = false;
    
    if(HasLattice(aTrack) == true){
        bIsUnderCoherentEffect = IsUnderCoherentEffect(aTrack);
        if(bIsUnderCoherentEffect){
            if(GetInfo(aTrack)->HasBeenUnderCoherentEffect() == 1){
                // if the particle is in channeling it gives the direction of the lattice to the particle momentum
                G4ThreeVector vPosition = ComputePositionInTheCrystal(aTrack.GetStep()->GetPostStepPoint(),aTrack);
                aParticleChange.ProposeMomentumDirection(GetXPhysicalLattice(aTrack)->GetLatticeDirection(vPosition).unit());
            }
            else if(GetInfo(aTrack)->HasBeenUnderCoherentEffect() == 2){
                // if the particle is in VR it gives a kick to the opposite site of the bending to the particle
                aParticleChange.ProposeMomentumDirection(ComputeVolumeReflectionOutgoingMomentum(aTrack));
            }
        }
    }
    else{
        // if the volume has no lattice it resets the density factors
        ResetDensity(aTrack);
    }
    
    if(bIsUnderCoherentEffect == false && (HasLattice(aTrack) == true || HasLatticeOnBoundary(aTrack) == true) ){
        // if has been under coherent effect but now it is not, the outgoing momentum is evaluated starting from the current position
        if(GetInfo(aTrack)->HasBeenUnderCoherentEffect() == 1){
            aParticleChange.ProposeMomentumDirection(ComputeChannelingOutgoingMomentum(aTrack));
            PrintChannelingPropertiesOnFile(aTrack);
        }
        else if(GetInfo(aTrack)->HasBeenUnderCoherentEffect() == 2){
            PrintChannelingPropertiesOnFile(aTrack);
        }
        
        // If is not under coherent effect sets coherent effect to zero and resets the density factors after the outgoing angle has been evaluated
            GetInfo(aTrack)->SetCoherentEffect(0);
        ResetDensity(aTrack);
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
    if(fLatticeManager->HasLattice(aTrack.GetStep()->GetPostStepPoint()->GetPhysicalVolume())){
        return fLatticeManager->GetXPhysicalLattice(aTrack.GetStep()->GetPostStepPoint()->GetPhysicalVolume());
    }
    else if(fLatticeManager->HasLattice(aTrack.GetStep()->GetPreStepPoint()->GetPhysicalVolume()) &&
            aTrack.GetStep()->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {
        return fLatticeManager->GetXPhysicalLattice(aTrack.GetStep()->GetPreStepPoint()->GetPhysicalVolume());
    }
    else{
        G4cout << "LATTICE NOT FOUND: ERROR" << G4endl;
        return NULL;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool ProcessChanneling::HasLattice(const G4Track& aTrack){
    if(fLatticeManager->HasLattice(aTrack.GetStep()->GetPostStepPoint()->GetPhysicalVolume())){
        return true;
    }
    else{
        return false;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool ProcessChanneling::HasLatticeOnBoundary(const G4Track& aTrack){
    if(fLatticeManager->HasLattice(aTrack.GetStep()->GetPreStepPoint()->GetPhysicalVolume()) &&
       aTrack.GetStep()->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {
        return true;
    }
    else{
        return false;
    }
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

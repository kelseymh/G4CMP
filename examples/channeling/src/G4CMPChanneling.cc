/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#include "G4CMPChanneling.hh"

#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4VParticleChange.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "G4RandomTools.hh"


#include "G4Box.hh"
#include "G4Tubs.hh"

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

G4CMPChanneling::G4CMPChanneling(const G4String& aName):G4VDiscreteProcess(aName){
    fLatticeManager = XLatticeManager3::GetXLatticeManager();
    
    G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
    G4cout<<"\n G4CMPChanneling::Constructor: Geometry surface tolerance is: " << kCarTolerance / mm << " mm"<<std::endl;
    if(verboseLevel>1) {
        G4cout << GetProcessName() << " is created "<< G4endl;
    }
    
    fFileName = "";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4CMPChanneling::~G4CMPChanneling(){
    fFileOut.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4CMPChanneling::G4CMPChanneling(G4CMPChanneling& right):G4VDiscreteProcess(right){
    ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CMPChanneling::SetFileName(const G4String& vFilename){
    
    if(fFileOut && vFilename && vFilename!=fFileName){
        fFileOut.close();
        
        fFileName = vFilename;
        
        G4String filename;
        fFileOut.open(filename = fFileName + "_CH.txt");
        fFileOut << "index,posin,angin,depth,pos,ang,angout,dens,tr_en,ndch,cheff" << std::endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4String G4CMPChanneling::GetFileName(){
    return fFileName;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CMPChanneling::PrintChannelingPropertiesOnFile(const G4Track& aTrack){
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

XVCrystalCharacteristic* G4CMPChanneling::GetPotential(){
    return fPotentialEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CMPChanneling::SetPotential(XVCrystalCharacteristic* vPotential){
    fPotentialEnergy = vPotential;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XVCrystalCharacteristic* G4CMPChanneling::GetElectricField(){
    return fElectricField;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CMPChanneling::SetElectricField(XVCrystalCharacteristic* vElectricField){
    fElectricField = vElectricField;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XCrystalIntegratedDensityHub* G4CMPChanneling::GetIntegratedDensity(){
    return fIntegratedDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CMPChanneling::SetIntegratedDensity(XCrystalIntegratedDensityHub* vIntegratedDensity){
    fIntegratedDensity = vIntegratedDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CMPChanneling::UpdateParameters(const G4Track& aTrack){
    
    if(fIntegratedDensity->HasBeenInitialized(GetXPhysicalLattice(aTrack)) == false){
        ComputeCrystalCharacteristic(aTrack);
        G4cout << "ChannelingProcess::UpdatePositionMomentumDensity::fIntegratedDensity->Initialized" << std::endl;
    }
    
    UpdatePosition(aTrack);
    UpdateMomentum(aTrack);
    UpdateDensity(aTrack);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CMPChanneling::UpdatePosition(const G4Track& aTrack){
    if(GetInfo(aTrack)->GetPositionChanneledInitial().x() == DBL_MAX || HasLatticeOnBoundary(aTrack)){
        G4double vXposition = G4UniformRand() * GetXPhysicalLattice(aTrack)->ComputeInterplanarPeriod();
        
        if(ParticleIsNegative(aTrack)){
            vXposition += ( GetXPhysicalLattice(aTrack)->ComputeInterplanarPeriod() * 0.5 );
        }
        
        GetInfo(aTrack)->SetPositionChanneled(G4ThreeVector(0.,0.,0.));
        
        GetInfo(aTrack)->SetPositionChanneledInitial(G4ThreeVector(vXposition,0.,0.));
    }
    else{
        G4double vPositionX = GetInfo(aTrack)->GetPositionChanneled().x();
        
        if(GetInfo(aTrack)->HasBeenUnderCoherentEffect() != 0 || HasLatticeOnBoundary(aTrack)){
            GetInfo(aTrack)->SetPositionChanneled(G4ThreeVector(0.,0.,0.));
        }
        else{
            vPositionX += (ComputePositionInTheCrystal(aTrack.GetStep()->GetPostStepPoint(),aTrack).r() - ComputePositionInTheCrystal(aTrack.GetStep()->GetPreStepPoint(),aTrack).r());
            GetInfo(aTrack)->SetPositionChanneled(G4ThreeVector(vPositionX,0.,0.));
        }        
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CMPChanneling::UpdateMomentum(const G4Track& aTrack){
    if(GetInfo(aTrack)->GetMomentumChanneledInitial().x() == DBL_MAX){
        // the first time it enter the crystal we take the momentum for the post step which is the only one in the crystal
        G4ThreeVector vMomentum = ComputeMomentum(aTrack,aTrack.GetStep()->GetPostStepPoint());
        
        GetInfo(aTrack)->SetMomentumChanneled(vMomentum);
        
        GetInfo(aTrack)->SetMomentumChanneledInitial(GetInfo(aTrack)->GetMomentumChanneled());
    }
    else{
        // we take the PREVIOUS step point to compare, otherwise the momentum is not computed correctly
       G4ThreeVector vMomentum = GetInfo(aTrack)->GetMomentumChanneled();
        
        vMomentum += ComputeMomentum(aTrack,aTrack.GetStep()->GetPreStepPoint());
        
        GetInfo(aTrack)->SetMomentumChanneled(vMomentum);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CMPChanneling::UpdateDensity(const G4Track& aTrack){
    
    
    G4double vTransverseEnergy = ComputeTransverseEnergy(aTrack).x();
    
    G4double vFactor = 1.;
        
    G4double vCharge = GetParticleDefinition(aTrack)->GetPDGCharge();
    G4double vNucleiDensity = fIntegratedDensity->GetIntegratedDensityNuclei(vTransverseEnergy,GetXPhysicalLattice(aTrack),vCharge);
    G4double vElectronDensity = fIntegratedDensity->GetIntegratedDensityElectron(vTransverseEnergy,GetXPhysicalLattice(aTrack),vCharge);
    
    GetInfo(aTrack)->SetNucleiDensity(vFactor * vNucleiDensity);
    GetInfo(aTrack)->SetElectronDensity(vFactor * vElectronDensity);
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CMPChanneling::ResetDensity(const G4Track& aTrack){
    GetInfo(aTrack)->SetNucleiDensity(1.);
    GetInfo(aTrack)->SetElectronDensity(1.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector G4CMPChanneling::ComputeChannelingOutgoingMomentum(const G4Track& aTrack){
    
    G4StepPoint* vStepPre = aTrack.GetStep()->GetPreStepPoint();
    G4StepPoint* vStepPost = aTrack.GetStep()->GetPostStepPoint();
    
    G4double vTotalEnergy = vStepPre->GetTotalEnergy();
    
    G4double vTransverseEnergyX = fabs(ComputeTransverseEnergy(aTrack).x() - ComputeCriticalEnergyMinimum(aTrack));
    G4double vTransverseEnergyZ = fabs(ComputeTransverseEnergy(aTrack).z() - ComputeCriticalEnergyMinimum(aTrack));

    G4double vChAngleX = pow(+ 2. * fabs(vTransverseEnergyX) / vTotalEnergy , 0.5);
    G4double vChAngleZ = pow(+ 2. * fabs(vTransverseEnergyZ) / vTotalEnergy , 0.5);
    
    G4double vPhi = 2. * ( G4UniformRand() - 0.5) * vChAngleX;
    G4double vTheta = 2. * ( G4UniformRand() - 0.5) * vChAngleZ;
        
    G4ThreeVector vNewMomentum = G4ThreeVector(0.,1.,0.).rotate(G4ThreeVector(0,0,1),- vPhi).rotate(G4ThreeVector(1,0,0),- vTheta);
    G4ThreeVector vPosition = ComputePositionInTheCrystal(vStepPost,aTrack);
    
    return GetXPhysicalLattice(aTrack)->ProjectMomentumVectorFromLatticeToWorld(vNewMomentum,vPosition);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector G4CMPChanneling::ComputeVolumeReflectionOutgoingMomentum(const G4Track& aTrack){
    
    G4StepPoint* vStep = aTrack.GetStep()->GetPostStepPoint();
    
    G4double vVrAngle = 0.;
    
    if(GetXPhysicalLattice(aTrack)->IsBent()) {
        G4double vRadiusX = GetXPhysicalLattice(aTrack)->GetCurvatureRadius().x();
        
        G4double vTotalEnergy = vStep->GetTotalEnergy();
        
        G4ThreeVector vInterplanarPeriod = G4ThreeVector(GetXPhysicalLattice(aTrack)->ComputeInterplanarPeriod(),0.,0.);
        
        G4double vTransverseEnergy = fabs(ComputeCriticalEnergyMaximum(aTrack) - ComputeCriticalEnergyMinimum(aTrack));
        
        G4double vCentrifugalForce = 0.;
        
        if(ParticleIsNegative(aTrack)){
            vCentrifugalForce = ComputeCentrifugalEnergy(aTrack, vInterplanarPeriod * 1.5).x();
        }
        else{
            vCentrifugalForce = ComputeCentrifugalEnergy(aTrack, vInterplanarPeriod).x();
        }
        
        vTransverseEnergy += (G4UniformRand() * fabs(vCentrifugalForce) );
        
        vVrAngle = - fabs(vRadiusX)/vRadiusX * pow(+ 2. * fabs(vTransverseEnergy) / vTotalEnergy , 0.5);
        
        if(ParticleIsNegative(aTrack)){
            vVrAngle *= 0.8; // = see PLB 681 (2009) 233
        }
        else{
            vVrAngle *= 1.4;   
        }
    }
    
    G4double vOmega = GetXPhysicalLattice(aTrack)->GetLatticeAngles().z();
    G4double vPhi = vVrAngle * cos(vOmega);
    G4double vTheta = vVrAngle * sin(vOmega);
    
    G4ThreeVector vNewMomentum = aTrack.GetMomentum().unit().rotate(G4ThreeVector(0,0,1), - vPhi).rotate(G4ThreeVector(0,1,0), -vTheta);
    
    return vNewMomentum.unit();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4ThreeVector G4CMPChanneling::ComputePositionInTheCrystal(G4StepPoint* vStep,const G4Track& aTrack){
    
    G4StepPoint* vStepVol = CheckStepPointLatticeForVolume(vStep,aTrack);
    G4StepPoint* vStepPos = CheckStepPointLatticeForPosition(vStep,aTrack);
    
    G4TouchableHistory* theTouchable = (G4TouchableHistory*)(vStepVol->GetTouchable());
    G4ThreeVector vWorldPos = vStepPos->GetPosition();
    G4ThreeVector vLocalPos = theTouchable->GetHistory()->GetTopTransform().TransformPoint(vWorldPos);
//    G4Tubs* vXtalSolid = (G4Tubs*) vStepVol->GetPhysicalVolume()->GetLogicalVolume()->GetSolid();
//    vLocalPos.rotateZ(-vXtalSolid->GetStartPhiAngle());
    G4Box* vXtalSolid = (G4Box*) vStepVol->GetPhysicalVolume()->GetLogicalVolume()->GetSolid();
    vLocalPos += G4ThreeVector(vXtalSolid->GetXHalfLength(),vXtalSolid->GetYHalfLength(),vXtalSolid->GetZHalfLength());
 
    return vLocalPos;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4StepPoint* G4CMPChanneling::CheckStepPointLatticeForVolume(G4StepPoint* vStep, const G4Track& aTrack){
    
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

G4StepPoint* G4CMPChanneling::CheckStepPointLatticeForPosition(G4StepPoint* vStep, const G4Track& aTrack){
    
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
        return vStepPre;
    }
    else{
        return vStep;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4CMPChanneling::IsUnderCoherentEffect(const G4Track& aTrack){
    //----------------------------------------
    // check if the particle momentum
    // transverse to the (h,k,l) plane
    // is small enough to permit channeling
    //----------------------------------------
        
    UpdateParameters(aTrack);
    
    G4double vEnergyMax = ComputeCriticalEnergyMaximum(aTrack);
    G4double vTransverseEnergy = ComputeTransverseEnergy(aTrack).x();
    
    if(GetXPhysicalLattice(aTrack)->IsBent() == false){
        
        if(vTransverseEnergy <= vEnergyMax){
            GetInfo(aTrack)->SetCoherentEffect(1);
            // the particle is in channeling
            return true;
        }
        else{
            // the particle is not under coherent effect
            GetInfo(aTrack)->SetCoherentEffect(0);
            return false;
        }
    }
    else{
        G4ThreeVector vPositionInTheCrystal = G4ThreeVector(GetInfo(aTrack)->GetPositionChanneled().x() + GetInfo(aTrack)->GetPositionChanneledInitial().x(),0.,0.);
        vTransverseEnergy += ComputeCentrifugalEnergy(aTrack,vPositionInTheCrystal).x();
        
        G4bool bNotBoundary = ParticleIsNotOnBoundary(aTrack);
        G4bool bCrossingPlane = ParticleIsCrossingPlane(aTrack);
        if(vTransverseEnergy <= vEnergyMax){
            // the particle is in channeling
            GetInfo(aTrack)->SetCoherentEffect(1);
            return true;
        }
        else if(bCrossingPlane == true &&
                bNotBoundary == true){
            // the particle is in volume reflection
            GetInfo(aTrack)->SetCoherentEffect(2);
            return true;
        }
        else{
            // the particle is not under coherent effect
            GetInfo(aTrack)->SetCoherentEffect(0);
            return false;
        }
    }
    return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4double G4CMPChanneling::GetChannelingMeanFreePath(const G4Track& aTrack){
    //----------------------------------------
    // return the channeling MFP
    //----------------------------------------
    
    G4double vMFP = ComputeOscillationPeriod(aTrack);
    
    if(GetInfo(aTrack)->GetNucleiDensity() < 1.){
        vMFP = ComputeOscillationPeriod(aTrack) / GetInfo(aTrack)->GetNucleiDensity();
    }
    
    return vMFP;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4CMPChanneling::GetMeanFreePath(const G4Track& aTrack,
                                            G4double, // previousStepSize
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

G4VParticleChange* G4CMPChanneling::PostStepDoIt(const G4Track& aTrack,
                                                   const G4Step&){
    
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
    
    if( (bIsUnderCoherentEffect == false && (HasLattice(aTrack) == true) ) || (HasLatticeOnBoundary(aTrack) == true) ) {
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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector G4CMPChanneling::ComputeTransverseEnergy(const G4Track& aTrack){
    //----------------------------------------
    // compute the particle transverse energy
    // in the crystal reference system
    //----------------------------------------
    
    G4ThreeVector vTransverseEnergy = ComputePotentialEnergy(aTrack) + ComputeKineticEnergy(aTrack);
    return vTransverseEnergy;
}

G4ThreeVector G4CMPChanneling::ComputeKineticEnergy(const G4Track& aTrack){
    //----------------------------------------
    // compute the particle kinetic energy
    // in the crystal reference system
    //----------------------------------------
    
    G4double vTotalEnergy = aTrack.GetStep()->GetPostStepPoint()->GetTotalEnergy();
    
    G4ThreeVector vMom = GetInfo(aTrack)->GetMomentumChanneled();
    
    G4ThreeVector vKineticEnergy = 0.5 * G4ThreeVector((vMom.x() * vMom.x()) / vTotalEnergy, 0., (vMom.z() * vMom.z())  / vTotalEnergy);
    return vKineticEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector G4CMPChanneling::ComputePotentialEnergy(const G4Track& aTrack){
    //----------------------------------------
    // compute the particle transverse energy
    // in the crystal reference system
    //----------------------------------------
    
    
    G4ThreeVector vPotentialEnergy = fPotentialEnergy->GetEC(GetInfo(aTrack)->GetPositionChanneledInitial(),GetXPhysicalLattice(aTrack));
    
    vPotentialEnergy *= GetParticleDefinition(aTrack)->GetPDGCharge();

    return vPotentialEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector G4CMPChanneling::ComputeCentrifugalEnergy(const G4Track& aTrack,G4ThreeVector vPosition){
    //----------------------------------------
    // compute the transverse energy variation
    // in the crystal reference system
    //----------------------------------------
    
    G4double vTotalEnergy = aTrack.GetStep()->GetPreStepPoint()->GetTotalEnergy();
    
    G4double vPositionX = vPosition.x();
    
    if(ParticleIsNegative(aTrack)){
        vPositionX -= GetXPhysicalLattice(aTrack)->ComputeInterplanarPeriod() * 0.5;
    }
    
    G4ThreeVector vEnergyVariation = G4ThreeVector(vTotalEnergy * vPositionX / GetXPhysicalLattice(aTrack)->GetCurvatureRadius().x(),0.,0.);
    
    return vEnergyVariation;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector G4CMPChanneling::ComputeMomentum(const G4Track& aTrack,G4StepPoint* vStep){
    //----------------------------------------
    // compute the particle momentum
    // in the crystal reference system
    //----------------------------------------
    
    G4ThreeVector vPosition = ComputePositionInTheCrystal(vStep,aTrack);
    
    G4ThreeVector vMomentum = GetXPhysicalLattice(aTrack)->ProjectMomentumVectorFromWorldToLattice(aTrack.GetMomentum(),vPosition);
    
    return vMomentum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4double G4CMPChanneling::ComputeCriticalEnergyMaximum(const G4Track& aTrack){
    //----------------------------------------
    // compute the critical energy
    // for channeling
    //----------------------------------------
    
    G4double vCriticalEnergy = 0.;
    
    if(ParticleIsNegative(aTrack)){
        vCriticalEnergy = - fPotentialEnergy->GetMinimum(GetXPhysicalLattice(aTrack));
    }
    else{
        vCriticalEnergy = + fPotentialEnergy->GetMaximum(GetXPhysicalLattice(aTrack));
    }
    
    vCriticalEnergy *= fabs(GetParticleDefinition(aTrack)->GetPDGCharge());
    
    return vCriticalEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4CMPChanneling::ComputeCriticalEnergyMinimum(const G4Track& aTrack){
    //----------------------------------------
    // compute the critical energy minimum
    // for channeling
    //----------------------------------------
    
    G4double vCriticalEnergy = 0.;
    
    if(ParticleIsNegative(aTrack)){
        vCriticalEnergy = - fPotentialEnergy->GetMaximum(GetXPhysicalLattice(aTrack));
    }
    else{
        vCriticalEnergy = + fPotentialEnergy->GetMinimum(GetXPhysicalLattice(aTrack));
    }
    
    vCriticalEnergy *= fabs(GetParticleDefinition(aTrack)->GetPDGCharge());
    
    return vCriticalEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4CMPChanneling::ComputeCriticalAngle(const G4Track& aTrack){
    //----------------------------------------
    // compute the critical angle
    // for chenneling
    //----------------------------------------
    
    G4double vTotalEnergy = aTrack.GetStep()->GetPostStepPoint()->GetTotalEnergy();
    G4double vCriticalAngle = pow( 2.0 * fabs( ( ComputeCriticalEnergyMaximum(aTrack) - ComputeCriticalEnergyMinimum(aTrack) ) / vTotalEnergy ) , 0.5);
    return vCriticalAngle;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4CMPChanneling::ComputeOscillationPeriod(const G4Track& aTrack){
    //----------------------------------------
    // compute the particle oscillation
    // period in the crystal channel
    //----------------------------------------
    
    G4double vInterplanarPeriod = GetXPhysicalLattice(aTrack)->ComputeInterplanarPeriod();
    G4double vOscillationPeriod = M_PI * vInterplanarPeriod / ComputeCriticalAngle(aTrack);
    return vOscillationPeriod;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4CMPChanneling::ComputeCriticalRadius(const G4Track& aTrack){
    //----------------------------------------
    // compute the critical radius
    // for channeling
    //----------------------------------------
    
    G4double vTotalEnergy = aTrack.GetStep()->GetPreStepPoint()->GetTotalEnergy();
    G4double vCriticalRadius = vTotalEnergy / fElectricField->GetMaximum(GetXPhysicalLattice(aTrack)) ;
    return vCriticalRadius;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector G4CMPChanneling::ComputePotentialWellCentre(const G4Track& aTrack){
    //----------------------------------------
    // compute the central point of
    // the potential well for channeling
    //----------------------------------------
    
    G4double vInterplanarPeriodHalf = 0.5 * GetXPhysicalLattice(aTrack)->ComputeInterplanarPeriod();
    
    G4double vCentreX = vInterplanarPeriodHalf;
    
    if(GetXPhysicalLattice(aTrack)->IsBent()){
        G4double vTotalEnergy = aTrack.GetStep()->GetPreStepPoint()->GetTotalEnergy();
        
        G4double vPotentialWellDepth = ComputeCriticalEnergyMaximum(aTrack) - (ComputeCriticalEnergyMinimum(aTrack));
        
        vCentreX *= (1. - 0.5 * vTotalEnergy / vPotentialWellDepth / GetXPhysicalLattice(aTrack)->GetCurvatureRadius().x() * vInterplanarPeriodHalf );
        
        if(ParticleIsNegative(aTrack)){
            vCentreX = (vInterplanarPeriodHalf * 2. - vCentreX);
        }
    }
    
    G4ThreeVector vCentre = G4ThreeVector(vCentreX,0.,0.);
    
    
    return vCentre;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4bool G4CMPChanneling::IsApplicable(const G4ParticleDefinition& aPD){
    return(aPD.GetPDGCharge() != 0.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CMPChanneling::BuildPhysicsTable(const G4ParticleDefinition&){
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XPhysicalLattice* G4CMPChanneling::GetXPhysicalLattice(const G4Track& aTrack){
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

G4bool G4CMPChanneling::HasLattice(const G4Track& aTrack){
    if(fLatticeManager->HasLattice(aTrack.GetStep()->GetPostStepPoint()->GetPhysicalVolume())){
        return true;
    }
    else{
        return false;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4CMPChanneling::HasLatticeOnBoundary(const G4Track& aTrack){
    if(fLatticeManager->HasLattice(aTrack.GetStep()->GetPreStepPoint()->GetPhysicalVolume()) &&
       aTrack.GetStep()->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {
        return true;
    }
    if(fLatticeManager->HasLattice(aTrack.GetStep()->GetPostStepPoint()->GetPhysicalVolume()) &&
       aTrack.GetStep()->GetPreStepPoint()->GetStepStatus() == fGeomBoundary) {
        return true;
    }
    else{
        return false;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4CMPChanneling::ParticleIsNegative(const G4Track& aTrack){
    if(GetParticleDefinition(aTrack)->GetPDGCharge() < 0.) {
        return true;
    }
    else{
        return false;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4CMPChanneling::ParticleIsCrossingPlane(const G4Track& aTrack){
    G4ThreeVector vPositionPre = ComputePositionInTheCrystal(aTrack.GetStep()->GetPreStepPoint(),aTrack);
    G4ThreeVector vMomentumPre = GetXPhysicalLattice(aTrack)->ProjectMomentumVectorFromWorldToLattice(aTrack.GetStep()->GetPreStepPoint()->GetMomentum(),vPositionPre);
    
    G4ThreeVector vPositionPost = ComputePositionInTheCrystal(aTrack.GetStep()->GetPostStepPoint(),aTrack);
    G4ThreeVector vMomentumPost = GetXPhysicalLattice(aTrack)->ProjectMomentumVectorFromWorldToLattice(aTrack.GetStep()->GetPostStepPoint()->GetMomentum(),vPositionPost);
    
    if(vMomentumPost.x()<0. &&
       vMomentumPre.x()>0.){
        return true;
    }
    if(vMomentumPost.x()>0. &&
       vMomentumPre.x()<0.){
        return true;
    }
    
    return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4CMPChanneling::ParticleIsNotOnBoundary(const G4Track& aTrack){
    if(aTrack.GetStep()->GetPostStepPoint()->GetStepStatus() != fGeomBoundary &&
       aTrack.GetStep()->GetPreStepPoint()->GetStepStatus() != fGeomBoundary){
        return true;
    }
    else{
        return false;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ChannelingParticleUserInfo* G4CMPChanneling::GetInfo(const G4Track& aTrack){
    return (ChannelingParticleUserInfo*) aTrack.GetUserInformation();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ParticleDefinition* G4CMPChanneling::GetParticleDefinition(const G4Track& aTrack){
    return const_cast<G4ParticleDefinition*>(aTrack.GetParticleDefinition());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CMPChanneling::ComputeCrystalCharacteristic(const G4Track& aTrack){
    fIntegratedDensity->SetXPhysicalLattice(GetXPhysicalLattice(aTrack));
    fIntegratedDensity->InitializeTables();

    fPotentialEnergy->InitializePhysicalLattice(GetXPhysicalLattice(aTrack));
    fElectricField->InitializePhysicalLattice(GetXPhysicalLattice(aTrack));

    if(fFileName != "") {
        PrintCrystalCharacteristicsOnFiles(aTrack);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CMPChanneling::PrintCrystalCharacteristicsOnFiles(const G4Track& aTrack){
    
    G4String filename;
    
    fIntegratedDensity->PrintOnFiles(filename=fFileName);
    fPotentialEnergy->PrintOnFile(filename=fFileName + "_pot.txt",GetXPhysicalLattice(aTrack),eV);
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

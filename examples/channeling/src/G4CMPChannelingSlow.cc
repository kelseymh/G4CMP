/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#include "G4CMPChannelingSlow.hh"

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

G4CMPChannelingSlow::G4CMPChannelingSlow(const G4String& aName):G4VDiscreteProcess(aName){
    fLatticeManager = XLatticeManager3::GetXLatticeManager();
    
    G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
    G4cout<<"\n G4CMPChannelingSlow::Constructor: Geometry surface tolerance is: " << kCarTolerance / mm << " mm"<<std::endl;
    if(verboseLevel>1) {
        G4cout << GetProcessName() << " is created "<< G4endl;
    }
    
    fFileName = "";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4CMPChannelingSlow::~G4CMPChannelingSlow(){
    fFileOut.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4CMPChannelingSlow::G4CMPChannelingSlow(G4CMPChannelingSlow& right):G4VDiscreteProcess(right){
    ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CMPChannelingSlow::SetFileName(const G4String& vFilename){
    
    if(fFileOut && vFilename && vFilename!=fFileName){
        fFileOut.close();
        
        fFileName = vFilename;
        
        G4String filename;
        fFileOut.open(filename = fFileName + "_CH.txt");
        fFileOut << "index,posin,angin,depth,pos,ang,angout,dens,tr_en,ndch,cheff" << std::endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4String G4CMPChannelingSlow::GetFileName(){
    return fFileName;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CMPChannelingSlow::PrintChannelingPropertiesOnFile(const G4Track& aTrack){
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

XVCrystalCharacteristic* G4CMPChannelingSlow::GetPotential(){
    return fPotentialEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CMPChannelingSlow::SetPotential(XVCrystalCharacteristic* vPotential){
    fPotentialEnergy = vPotential;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XVCrystalCharacteristic* G4CMPChannelingSlow::GetElectricField(){
    return fElectricField;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CMPChannelingSlow::SetElectricField(XVCrystalCharacteristic* vElectricField){
    fElectricField = vElectricField;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XCrystalIntegratedDensityHub* G4CMPChannelingSlow::GetIntegratedDensity(){
    return fIntegratedDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CMPChannelingSlow::SetIntegratedDensity(XCrystalIntegratedDensityHub* vIntegratedDensity){
    fIntegratedDensity = vIntegratedDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CMPChannelingSlow::UpdateParameters(const G4Track& aTrack){
    
    if(fIntegratedDensity->HasBeenInitialized(GetXPhysicalLattice(aTrack)) == false){
        ComputeCrystalCharacteristic(aTrack);
        G4cout << "ChannelingProcess::UpdateParameters::fIntegratedDensity->Initialized" << std::endl;
    }
    
    UpdatePosition(aTrack);
    UpdateMomentum(aTrack);
    
    G4double vStepLengthTotal = aTrack.GetStepLength();
    G4double vStepLengthMin = 1. * CLHEP::angstrom;
    G4double vStepLengthMax = 0.01 * CLHEP::micrometer;
    G4double vTransverseVariationMax = 1.E-2 * CLHEP::angstrom;
    
    G4double vStepLength = 0.;
    G4double vStepLengthHalf = 0.;

    G4double vCharge = GetParticleDefinition(aTrack)->GetPDGCharge();
    G4double vDensityNuclei = 0.;
    G4double vDensityElectron = 0.;
    G4int vNumberOfSteps = 0;
    
    G4double vPositionX = GetInfo(aTrack)->GetPositionChanneled().x();
    G4double vPositionY = 0.;
    G4double vPositionZ = GetInfo(aTrack)->GetPositionChanneled().z();
    G4double vMomentumX = GetInfo(aTrack)->GetMomentumChanneled().x();
    G4double vMomentumY = GetInfo(aTrack)->GetMomentumChanneled().y();
    G4double vMomentumZ = GetInfo(aTrack)->GetMomentumChanneled().z();
    G4double vPositionHalfX = GetInfo(aTrack)->GetPositionChanneled().x();
    G4double vMomentumHalfX = GetInfo(aTrack)->GetMomentumChanneled().x();
    
    G4double vRadiusX = GetXPhysicalLattice(aTrack)->GetCurvatureRadius().x();

    G4double vInterplanarPeriod = GetXPhysicalLattice(aTrack)->ComputeInterplanarPeriod();

    vPositionX = 0.7 * CLHEP::angstrom;
    G4bool bExit = false;
    do{
        
        if(vMomentumX != 0.0){
            vStepLength = fabs(vTransverseVariationMax / vMomentumX * vMomentumY);
            if( vStepLength <= vStepLengthMin ) vStepLength = vStepLengthMin;
            if( vStepLength >= vStepLengthMax ) vStepLength = vStepLengthMax;
        }
        else{
            vStepLength = vStepLengthMin;
        }
        
        if(vPositionY>vStepLengthTotal){
            vStepLength = vStepLengthTotal - vPositionY;
            bExit = true;
        }
        
        vStepLengthHalf = vStepLength * 0.5;
        
        vPositionHalfX += vMomentumX / vMomentumY * vStepLengthHalf;
        vMomentumHalfX += vCharge * fElectricField->GetEC(G4ThreeVector(vPositionX,0.,0.),GetXPhysicalLattice(aTrack)).x() * vStepLengthHalf;
        
        if(vRadiusX!=0.){
            vMomentumHalfX -= ( vMomentumY * vStepLengthHalf / vRadiusX );
        }
        
        vPositionX += vMomentumHalfX / vMomentumY * vStepLength;
        vMomentumX += vCharge * fElectricField->GetEC(G4ThreeVector(vPositionHalfX,0.,0.),GetXPhysicalLattice(aTrack)).x() * vStepLength;
        
        if(vRadiusX!=0.){
            vMomentumX -= (vMomentumY * vStepLength / vRadiusX );
        }
        vPositionZ += vMomentumZ / vMomentumY * vStepLength;
        
        vPositionY += vStepLength;
        
        vDensityNuclei += fIntegratedDensity->GetDensityNuclei()->GetEC(G4ThreeVector(vPositionX,0.,0.),GetXPhysicalLattice(aTrack)).x();
        vDensityElectron += fIntegratedDensity->GetDensityElectron()->GetEC(G4ThreeVector(vPositionX,0.,0.),GetXPhysicalLattice(aTrack)).x();

        vNumberOfSteps++;
        
    } while(!bExit);
    
    GetInfo(aTrack)->SetMomentumChanneled(G4ThreeVector(vMomentumX,vMomentumY,vMomentumZ));
    GetInfo(aTrack)->SetPositionChanneled(G4ThreeVector(vPositionX,vPositionY,vPositionZ));
    GetInfo(aTrack)->SetNucleiDensity(vDensityNuclei * vInterplanarPeriod / vNumberOfSteps);
    GetInfo(aTrack)->SetElectronDensity(vDensityElectron * vInterplanarPeriod / vNumberOfSteps);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CMPChannelingSlow::UpdatePosition(const G4Track& aTrack){
    if(GetInfo(aTrack)->GetPositionChanneledInitial().x() == DBL_MAX || HasLatticeOnBoundary(aTrack)){
        G4double vXposition = G4UniformRand() * GetXPhysicalLattice(aTrack)->ComputeInterplanarPeriod();
        
        GetInfo(aTrack)->SetPositionChanneled(G4ThreeVector(0.,0.,0.));
        
        GetInfo(aTrack)->SetPositionChanneledInitial(G4ThreeVector(vXposition,0.,0.));
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CMPChannelingSlow::UpdateMomentum(const G4Track& aTrack){
    if(GetInfo(aTrack)->GetMomentumChanneledInitial().x() == DBL_MAX){
        // the first time it enter the crystal we take the momentum for the post step which is the only one in the crystal
        G4ThreeVector vMomentum = ComputeMomentum(aTrack,aTrack.GetStep()->GetPostStepPoint());
        
        GetInfo(aTrack)->SetMomentumChanneled(vMomentum);
        
        GetInfo(aTrack)->SetMomentumChanneledInitial(GetInfo(aTrack)->GetMomentumChanneled());
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CMPChannelingSlow::ResetDensity(const G4Track& aTrack){
    GetInfo(aTrack)->SetNucleiDensity(1.);
    GetInfo(aTrack)->SetElectronDensity(1.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4ThreeVector G4CMPChannelingSlow::ComputePositionInTheCrystal(G4StepPoint* vStep,const G4Track& aTrack){
    
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

G4StepPoint* G4CMPChannelingSlow::CheckStepPointLatticeForVolume(G4StepPoint* vStep, const G4Track& aTrack){
    
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

G4StepPoint* G4CMPChannelingSlow::CheckStepPointLatticeForPosition(G4StepPoint* vStep, const G4Track& aTrack){
    
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


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4CMPChannelingSlow::GetChannelingMeanFreePath(const G4Track& aTrack){
    //----------------------------------------
    // return the channeling MFP
    //----------------------------------------
    
    G4double vMFP = ComputeOscillationPeriod(aTrack);
        
    return vMFP;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4CMPChannelingSlow::GetMeanFreePath(const G4Track& aTrack,
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

G4VParticleChange* G4CMPChannelingSlow::PostStepDoIt(const G4Track& aTrack,
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
    
    if(HasLattice(aTrack) == true && !HasLatticeOnBoundary(aTrack)){
        UpdateParameters(aTrack);
        
        G4ThreeVector vPosition = ComputePositionInTheCrystal(aTrack.GetStep()->GetPostStepPoint(),aTrack);

        G4ThreeVector vMomentum = GetXPhysicalLattice(aTrack)->ProjectMomentumVectorFromLatticeToWorld(GetInfo(aTrack)->GetMomentumChanneled(),vPosition);

        aParticleChange.ProposeMomentumDirection(GetInfo(aTrack)->GetMomentumChanneled().unit());
        //G4cout << GetInfo(aTrack)->GetMomentumChanneled().x()/GetInfo(aTrack)->GetMomentumChanneled().y() << G4endl;
    }
    else{
        // if the volume has no lattice it resets the density factors
        ResetDensity(aTrack);
    }
    
    return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector G4CMPChannelingSlow::ComputeTransverseEnergy(const G4Track& aTrack){
    //----------------------------------------
    // compute the particle transverse energy
    // in the crystal reference system
    //----------------------------------------
    
    G4ThreeVector vTransverseEnergy = ComputePotentialEnergy(aTrack) + ComputeKineticEnergy(aTrack);
    return vTransverseEnergy;
}

G4ThreeVector G4CMPChannelingSlow::ComputeKineticEnergy(const G4Track& aTrack){
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

G4ThreeVector G4CMPChannelingSlow::ComputePotentialEnergy(const G4Track& aTrack){
    //----------------------------------------
    // compute the particle transverse energy
    // in the crystal reference system
    //----------------------------------------
    
    
    G4ThreeVector vPotentialEnergy = fPotentialEnergy->GetEC(GetInfo(aTrack)->GetPositionChanneledInitial(),GetXPhysicalLattice(aTrack));
    
    vPotentialEnergy *= GetParticleDefinition(aTrack)->GetPDGCharge();
    
    return vPotentialEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector G4CMPChannelingSlow::ComputeCentrifugalEnergy(const G4Track& aTrack,G4ThreeVector vPosition){
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

G4ThreeVector G4CMPChannelingSlow::ComputeMomentum(const G4Track& aTrack,G4StepPoint* vStep){
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


G4double G4CMPChannelingSlow::ComputeCriticalEnergyMaximum(const G4Track& aTrack){
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

G4double G4CMPChannelingSlow::ComputeCriticalEnergyMinimum(const G4Track& aTrack){
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

G4double G4CMPChannelingSlow::ComputeCriticalAngle(const G4Track& aTrack){
    //----------------------------------------
    // compute the critical angle
    // for chenneling
    //----------------------------------------
    
    G4double vTotalEnergy = aTrack.GetStep()->GetPostStepPoint()->GetTotalEnergy();
    G4double vCriticalAngle = pow( 2.0 * fabs( ( ComputeCriticalEnergyMaximum(aTrack) - ComputeCriticalEnergyMinimum(aTrack) ) / vTotalEnergy ) , 0.5);
    return vCriticalAngle;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4CMPChannelingSlow::ComputeOscillationPeriod(const G4Track& aTrack){
    //----------------------------------------
    // compute the particle oscillation
    // period in the crystal channel
    //----------------------------------------
    
    G4double vInterplanarPeriod = GetXPhysicalLattice(aTrack)->ComputeInterplanarPeriod();
    G4double vOscillationPeriod = M_PI * vInterplanarPeriod / ComputeCriticalAngle(aTrack);
    return vOscillationPeriod;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4CMPChannelingSlow::ComputeCriticalRadius(const G4Track& aTrack){
    //----------------------------------------
    // compute the critical radius
    // for channeling
    //----------------------------------------
    
    G4double vTotalEnergy = aTrack.GetStep()->GetPreStepPoint()->GetTotalEnergy();
    G4double vCriticalRadius = vTotalEnergy / fElectricField->GetMaximum(GetXPhysicalLattice(aTrack)) ;
    return vCriticalRadius;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector G4CMPChannelingSlow::ComputePotentialWellCentre(const G4Track& aTrack){
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


G4bool G4CMPChannelingSlow::IsApplicable(const G4ParticleDefinition& aPD){
    return(aPD.GetPDGCharge() != 0.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CMPChannelingSlow::BuildPhysicsTable(const G4ParticleDefinition&){
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XPhysicalLattice* G4CMPChannelingSlow::GetXPhysicalLattice(const G4Track& aTrack){
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

G4bool G4CMPChannelingSlow::HasLattice(const G4Track& aTrack){
    if(fLatticeManager->HasLattice(aTrack.GetStep()->GetPostStepPoint()->GetPhysicalVolume())){
        return true;
    }
    else{
        return false;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4CMPChannelingSlow::HasLatticeOnBoundary(const G4Track& aTrack){
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

G4bool G4CMPChannelingSlow::ParticleIsNegative(const G4Track& aTrack){
    if(GetParticleDefinition(aTrack)->GetPDGCharge() < 0.) {
        return true;
    }
    else{
        return false;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4CMPChannelingSlow::ParticleIsCrossingPlane(const G4Track& aTrack){
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

G4bool G4CMPChannelingSlow::ParticleIsNotOnBoundary(const G4Track& aTrack){
    if(aTrack.GetStep()->GetPostStepPoint()->GetStepStatus() != fGeomBoundary &&
       aTrack.GetStep()->GetPreStepPoint()->GetStepStatus() != fGeomBoundary){
        return true;
    }
    else{
        return false;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ChannelingParticleUserInfo* G4CMPChannelingSlow::GetInfo(const G4Track& aTrack){
    return (ChannelingParticleUserInfo*) aTrack.GetUserInformation();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ParticleDefinition* G4CMPChannelingSlow::GetParticleDefinition(const G4Track& aTrack){
    return const_cast<G4ParticleDefinition*>(aTrack.GetParticleDefinition());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CMPChannelingSlow::ComputeCrystalCharacteristic(const G4Track& aTrack){
    fIntegratedDensity->SetXPhysicalLattice(GetXPhysicalLattice(aTrack));
    fIntegratedDensity->InitializeTables();
    
    fPotentialEnergy->InitializePhysicalLattice(GetXPhysicalLattice(aTrack));
    fElectricField->InitializePhysicalLattice(GetXPhysicalLattice(aTrack));
    
    if(fFileName != "") {
        PrintCrystalCharacteristicsOnFiles(aTrack);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CMPChannelingSlow::PrintCrystalCharacteristicsOnFiles(const G4Track& aTrack){
    
    G4String filename;
    
    fIntegratedDensity->PrintOnFiles(filename=fFileName);
    fPotentialEnergy->PrintOnFile(filename=fFileName + "_pot.txt",GetXPhysicalLattice(aTrack),eV);
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

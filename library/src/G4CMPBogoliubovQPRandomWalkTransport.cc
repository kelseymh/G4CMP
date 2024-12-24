/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPBogoliubovQPRandomWalkTransport.hh
/// \brief Implementation of the G4CMPBogoliubovQPRandomWalkTransport class

#include "G4CMPBogoliubovQPRandomWalkTransport.hh"
#include "G4CMPSCUtils.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPVProcess.hh"
#include "G4CMPGeometryUtils.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4CMPProcessUtils.hh"
#include "G4CMPTrackUtils.hh"
#include "G4CMPVProcess.hh"
#include "G4TransportationManager.hh"
#include "G4SafetyHelper.hh"
#include "Randomize.hh"

#include <iostream>
#include <cmath>
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4CMPBogoliubovQPRandomWalkTransport::G4CMPBogoliubovQPRandomWalkTransport(const G4String& name , G4CMPProcessSubType fBogoliubovQPRandomWalkTransport)
  : G4VContinuousDiscreteProcess(name,fPhonon),G4CMPBoundaryUtils(this),
  fNewPosition(0.,0.,0.),
  fNewDirection(1.,0.,0.)
{
  verboseLevel = G4CMPConfigManager::GetVerboseLevel();
  SetProcessSubType(fBogoliubovQPRandomWalkTransport);
    
  //setting default time step to 10 us should be overwritten by SC lattice params
  fTimeStep = 10*us;

  //Default step values
  fStepX =  0.0;
  fStepY =  0.0;
  fPathLength =  0.0;
  fDiffConst =  0.0;

  //fSafetyHelper is initialized in AlongStepGPIL
  fSafetyHelper=nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4CMPBogoliubovQPRandomWalkTransport::~G4CMPBogoliubovQPRandomWalkTransport()
{
  if(verboseLevel>2) {
    G4cout << "G4CMPBogoliubovQPRandomWalkTransport destruct " << GetProcessName()
          << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CMPBogoliubovQPRandomWalkTransport::StartTracking(G4Track* track)
{
    G4VProcess::StartTracking(track);    // Apply base class actions
    LoadDataForTrack(track);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4CMPBogoliubovQPRandomWalkTransport::AlongStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double,
                             G4double currentMinimalStep,
                             G4double&,
                             G4GPILSelection* selection)
{
    if(!fSafetyHelper) {
        G4TransportationManager* transportMgr ;

        transportMgr = G4TransportationManager::GetTransportationManager() ;

        fSafetyHelper = transportMgr->GetSafetyHelper();
        
        fSafetyHelper->InitialiseHelper();
    }
    
    // get Step limit proposed by the post-step processes
    *selection = NotCandidateForSelection;
    fPathLength = currentMinimalStep;
    
    // Do not need to check the material table since this is loaded in via
    // the logical and physical lattice
    // If lattice is a superconductor then fDn != 0 and fGapEnergy != 0
    G4double energy = track.GetKineticEnergy();
    
    //Get velocity of the current track
    G4double velocity = track.GetVelocity();
    
    //Initialize the lattice manager
    G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();
    
    //Seems like it is probable for the pre and post step of the track to be same the volume if they are are
    //at the same position
    //This is a hacky way of checking if there is a lattice in the momentum direction passed to the AlongStepGPIL
    G4ThreeVector trackPos = track.GetPosition();
    G4ThreeVector pDir = track.GetMomentumDirection();
    G4ThreeVector trackPos_eps = trackPos+pDir*1*nm;
    G4ThreeVector trackPos_minus_eps = trackPos-pDir*1*nm;
    G4VPhysicalVolume* currentVol = G4CMP::GetVolumeAtPoint(trackPos);
    G4VPhysicalVolume* currentVol_eps = G4CMP::GetVolumeAtPoint(trackPos_eps);
    G4VPhysicalVolume* currentVol_minus_eps = G4CMP::GetVolumeAtPoint(trackPos_minus_eps);
    G4bool positionPlusEpsLattice =LM->HasLattice(currentVol_eps);
    G4bool positionMinusEpsLattice =LM->HasLattice(currentVol_minus_eps);
    
    if(verboseLevel>2) {
        G4cout << "current position has lattice status "<<LM->HasLattice(currentVol)<<G4endl;
        G4cout << "current position plus eps has lattice status "<<positionPlusEpsLattice<<G4endl;
        G4cout << "current position minus eps has lattice status "<<positionMinusEpsLattice<<G4endl;
    }
    
    //Check if the current step is on the boundary
    //it is likely that the below boolean can just check if the current steppoint has fGeomBoundary status
    fTrackOnBoundary = currentVol_eps != currentVol_minus_eps;
    
    G4ThreeVector surfaceNorm;
    
    if (fTrackOnBoundary){
        surfaceNorm=G4CMP::GetSurfaceNormal(*track.GetStep());
        if(verboseLevel>2) {
            G4cout << "G4CMPBogoliubovQPRandomWalkTransport: inside AlongStepGPIL : fTrackOnBoundary is true " <<"surfaceNorm : "<<surfaceNorm << G4endl;
            G4cout << "G4CMPBogoliubovQPRandomWalkTransport: inside AlongStepGPIL : Step status fGeomBoundary? : " << (track.GetStep()->GetPreStepPoint()->GetStepStatus()==fGeomBoundary) << G4endl;
            G4cout << "G4CMPBogoliubovQPRandomWalkTransport: inside AlongStepGPIL " <<"momentum : "<<pDir << G4endl;
            G4cout << "G4CMPBogoliubovQPRandomWalkTransport: inside AlongStepGPIL " <<"momentum dot surfaceNorm : "<<pDir.dot(surfaceNorm) << G4endl;
        }
    }
    
    G4double proposedTimeStep = 0.0;
    
    G4LatticePhysical* theLat;
    
    if (LM->HasLattice(currentVol)){
        //Current position volume has a lattice attached
        theLat = LM->GetLattice(currentVol);
        fGapEnergy = theLat->GetSCDelta0();
        fTau0_qp = theLat->GetSCTau0qp();
        fDn = theLat->GetSCDn();
        proposedTimeStep = 0.5*fTau0_qp;
    } else{
        //The current volume doesnt have a lattice attached return max double and make process inactive
        if(verboseLevel>2) G4cout << "G4CMPBogoliubovQPRandomWalkTransport: inside AlongStepGPIL " <<"currentVol does not have a lattice therefore making along step inactive " << G4endl;
        fStepX = 0.0;
        fStepY = 0.0;
        isActive = false;
        return DBL_MAX;
    }
    
    if ((energy>=fGapEnergy) && (fDn)>0){
        isActive = true;
        *selection = CandidateForSelection;
        
        //Calculate the energy dependent diffusion constant
        G4double E_ratio = fGapEnergy/energy;
//        G4cout << "Energy QP : " << energy <<G4endl;
//        G4cout << "Material gap : " << fGapEnergy <<G4endl;
//        G4cout << "E_ratio : " << E_ratio <<G4endl;
//        G4cout << "Normal State Diffusion Constant : " << fDn <<G4endl;
        G4double E_ratio2 = pow(E_ratio,2.0);
        fDiffConst = fDn*sqrt(1-E_ratio2);
//        G4cout << "QP effective Diffusion constant : " << fDiffConst <<G4endl;

        //Want to update update time step to be half of the phonon
        //scattering lifetime - estimate from Eq.5 in Saving paper Martinis 2021
        G4double E_ratio_diff = energy/fGapEnergy - 1;
        proposedTimeStep = 0.5*fTau0_qp/(1.8*pow(E_ratio_diff,3));
        
        //Time step is set by the current minimum process else it is half the scattering rate
//        G4cout << "Time step from pathlength and velocity : " << (fPathLength/velocity)*ns <<G4endl;
//        G4cout << "Track local time : " << track.GetLocalTime()*ns <<G4endl;
//        G4cout << "Track global time : " << track.GetLocalTime()*ns <<G4endl;
        fTimeStep = std::min(proposedTimeStep, fPathLength/velocity);
        
//        G4cout << "Taken time step : " << fTimeStep*ns <<G4endl;
        
        //Now we can calculate the std of the guassian for the RW
        G4double sigma = sqrt(2.0*fDiffConst*fTimeStep);
        
        G4RandGauss* gauss_dist = new G4RandGauss(G4Random::getTheEngine(),0.0,sigma);
        
//        fStepX = gauss_dist->shoot();
//        fStepY = gauss_dist->shoot();
        
        fStepX = gauss_dist->fire();
        fStepY = gauss_dist->fire();

        
        fPathLength = sqrt((fStepX*fStepX+fStepY*fStepY));
        //If the walker is on the boundary we want to make sure that
        //the proposed step agrees with the momentum direction proposed by the
        //QP boundary process class
        
        if (fTrackOnBoundary){
            //get the surface normal for picking the next postion of the random walk
            G4double pDotSurfaceNorm = pDir.dot(surfaceNorm);
            G4bool pDotSurfaceNormPositive = pDotSurfaceNorm>0;
            
            if(verboseLevel>2){
                G4cout << "G4CMPBogoliubovQPRandomWalkTransport: inside AlongStepGPIL " <<" is pDotSurfaceNorm positive? : "<<pDotSurfaceNormPositive<< G4endl;
                G4cout << "G4CMPBogoliubovQPRandomWalkTransport: inside AlongStepGPIL " <<"pDotSurfaceNorm : "<<pDotSurfaceNorm<< G4endl;
            }
            
            if (positionPlusEpsLattice){
                //So the position plus small delta along direction of momentum has a lattice
                // this would be the case if the momentum direction was reflected/transmitted by the boundary class
                if (pDotSurfaceNormPositive){
                    //If the momentum direction projected on the outward surface normal is positive then the particle has been transmitted according to the boundary class -> we must constrain the
                    //suggested steps wrt to this boundary
                    G4double phi = pi*G4UniformRand()-pi/2;
                    if(verboseLevel>2)G4cout << "G4CMPBogoliubovQPRandomWalkTransport: inside AlongStepGPIL " <<"phi for positive dot product : "<<phi<< G4endl;
                    fStepX = fPathLength*(surfaceNorm.x()*std::cos(phi)-surfaceNorm.y()*std::sin(phi));
                    fStepY = fPathLength*(surfaceNorm.x()*std::sin(phi)+surfaceNorm.y()*std::cos(phi));
                }
                else{
                    //If the momentum direction is projected on the outward surface normal is negative then the particles has been reflected by the boundary class ->
                    //We must constrain the suggested steps to be in this direction wrt the boundart
                    G4double phi = pi*G4UniformRand()-pi/2;
                    if(verboseLevel>2)G4cout << "G4CMPBogoliubovQPRandomWalkTransport: inside AlongStepGPIL " <<"phi for negative dot product : "<<phi<< G4endl;
                    fStepX = fPathLength*(-surfaceNorm.x()*std::cos(phi)+surfaceNorm.y()*std::sin(phi));
                    fStepY = fPathLength*(-surfaceNorm.x()*std::sin(phi)-surfaceNorm.y()*std::cos(phi));
                }
            }else if (positionMinusEpsLattice){
                //This case is handling the case where only the negative momentum direction has a lattice
                //We then want to treat the negative momentum direction as the true momentum direction
                if (pDotSurfaceNormPositive){
                    //If the momentum direction projected on the outward surface normal is positive then the particle has been transmitted according to the boundary class -> we must constrain the
                    //suggested steps wrt to this boundary
                    G4double phi = pi*G4UniformRand()-pi/2;
                    fStepX = fPathLength*(-surfaceNorm.x()*std::cos(phi)+surfaceNorm.y()*std::sin(phi));
                    fStepY = fPathLength*(-surfaceNorm.x()*std::sin(phi)-surfaceNorm.y()*std::cos(phi));
                }
                else{
                    //If the momentum direction is projected on the outward surface normal is negative then the particles has been reflected by the boundary class ->
                    //We must constrain the suggested steps to be in this direction wrt the bound
                    G4double phi = pi*G4UniformRand()-pi/2;
                    fStepX = fPathLength*(surfaceNorm.x()*std::cos(phi)-surfaceNorm.y()*std::sin(phi));
                    fStepY = fPathLength*(surfaceNorm.x()*std::sin(phi)+surfaceNorm.y()*std::cos(phi));
                }
            }
        }
        
        if(verboseLevel>1)G4cout << "G4CMPBogoliubovQPRandomWalkTransport: inside AlongStepGPIL " <<"x step : "<<fStepX <<"y step : "<<fStepY<< G4endl;
        
        return fPathLength;
    }
    else{
        fStepX = 0.0;
        fStepY = 0.0;
        isActive = false;
        return DBL_MAX;
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4CMPBogoliubovQPRandomWalkTransport::PostStepGetPhysicalInteractionLength(
              const G4Track&, G4double, G4ForceCondition* condition)
{
  *condition = NotForced;
  return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4CMPBogoliubovQPRandomWalkTransport::AlongStepDoIt(const G4Track& track, const G4Step& step) {
    
  fParticleChange.ProposeMomentumDirection(
    step.GetPostStepPoint()->GetMomentumDirection());
  fNewPosition = step.GetPostStepPoint()->GetPosition();
  fParticleChange.ProposePosition(fNewPosition);
  G4double velocity = step.GetPostStepPoint()->GetVelocity();
  fParticleChange.ProposeVelocity(velocity);
  fPositionChanged = false;

  G4double stepLength = step.GetStepLength();
    
  // Check if the particle met conditions to do random walk from GPIL command
  if(!isActive) {
    fPathLength = stepLength;
    G4cout << "Step is not active" << G4endl;
    G4cout << "Particle path length: " << fPathLength <<G4endl;
    G4cout << "Particle velocity: " << velocity <<G4endl;
    G4cout << "Time Step : " << fPathLength/velocity <<G4endl;
  } else {
    // Particle did meet conditions to undergo RW
    // Update relevant particle changes
    //'Velocity' of particle is set by displacement and simulation time step
    velocity = fPathLength/fTimeStep;
    fParticleChange.ProposeVelocity(velocity);
      
    //'MomentumDirection' of particle is set by the random walk gaussian steps
    fNewDirection = G4ThreeVector(fStepX/fPathLength,fStepY/fPathLength,0.0);
    fOldPosition = step.GetPreStepPoint()->GetPosition();
    fNewPosition = fOldPosition+fPathLength*fNewDirection;
      
    G4double dist = 0.0;
    G4double safety = fSafetyHelper->ComputeSafety(fOldPosition);
      
    if(fSafetyHelper->RecheckDistanceToCurrentBoundary(fOldPosition, fNewDirection, fPathLength, &dist, &safety)){
        if (fPathLength <=dist){
            if ((fPathLength==dist)&&fTrackOnBoundary){
                if(verboseLevel>2) G4cout << "G4CMPBogoliubovQPRandomWalkTransport: inside AlongStepDoIt \t fPathLength == dist : "<<dist << G4endl;
                fPositionChanged = false;
                fNewPosition=fOldPosition;
                fParticleChange.ProposePosition(fNewPosition);
            }else{
                if(verboseLevel>2){
                    G4cout << "G4CMPBogoliubovQPRandomWalkTransport: inside AlongStepDoIt \t dist : "<<dist << G4endl;
                    G4cout << "\t fPathLength : "<<fPathLength<<G4endl;
                    G4cout << "\t new position from GPIL : "<< fNewPosition << G4endl;
                    G4cout << "\t proposed new position : "<< fNewPosition<<G4endl;
                    G4cout << "\t new position : "<< fNewPosition << G4endl;
                }
                fSafetyHelper->ReLocateWithinVolume(fNewPosition);
                fParticleChange.ProposeMomentumDirection(fNewDirection);
                fPositionChanged = true;
            }
        }else{
            if(verboseLevel>2){
                G4cout << "G4CMPBogoliubovQPRandomWalkTransport: inside AlongStepDoIt  \t dist : "<<dist << G4endl;
                G4cout << "\t new position from GPIL : "<< fNewPosition << G4endl;
                G4cout << "\t proposed new position : "<< fOldPosition+(dist)*fNewDirection<<G4endl;
            }
            fNewPosition = fOldPosition+(dist)*fNewDirection;
            if(verboseLevel>2) G4cout << "\t new position : "<< fNewPosition << G4endl;
            fParticleChange.ProposeMomentumDirection(fNewDirection);
            fPositionChanged = true;
            fPathLength=dist;
        }
        
    }
      
   if(verboseLevel>2){
       G4cout << "G4CMPBogoliubovQPRandomWalkTransport: inside AlongStepDoIt  \t Old position : "<< fOldPosition <<G4endl;
       G4cout << "\t New position : "<< fNewPosition <<G4endl;
       G4cout << "\t New safety : "<< safety <<G4endl;
    }
    //the proposed local time is not updated by the UpdateAlongStep method in the particle change class
    // - > this is likely not needed and updated by the G4Transportation class
    fParticleChange.ProposeLocalTime(fPathLength/velocity);
  }
  if(fPositionChanged){
      fParticleChange.ProposePosition(fNewPosition);
  }
    
  fParticleChange.ProposeTrueStepLength(fPathLength);
    
//  G4cout << "Particle path length: " << fPathLength <<G4endl;
//  G4cout << "Particle velocity: " << velocity <<G4endl;
//  G4cout << "Time Step : " << fPathLength/velocity <<G4endl;
    
  return &fParticleChange;
}
    
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* 
G4CMPBogoliubovQPRandomWalkTransport::PostStepDoIt(const G4Track& track, const G4Step&)
{
  fParticleChange.Initialize(track);
  //ReloadDataForTrack(&track); //REL commenting this out because I no longer use this function and it will be deleted from the final branch
  return &fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4CMPBogoliubovQPRandomWalkTransport::GetContinuousStepLimit(
                                       const G4Track& track,
                                       G4double previousStepSize,
                                       G4double currentMinimalStep,
                                       G4double& currentSafety)
{
  G4GPILSelection selection = NotCandidateForSelection;
  G4double x = AlongStepGetPhysicalInteractionLength(track,previousStepSize,
                                                     currentMinimalStep,
                                                     currentSafety,
                                                     &selection);
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4CMPBogoliubovQPRandomWalkTransport::ContinuousStepLimit(
                                       const G4Track& track,
                                       G4double previousStepSize,
                                       G4double currentMinimalStep,
                                       G4double& currentSafety)
{
  return GetContinuousStepLimit(track,previousStepSize,currentMinimalStep,
                                currentSafety);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4CMPBogoliubovQPRandomWalkTransport::GetMeanFreePath(
              const G4Track&, G4double, G4ForceCondition* condition)
{
//  if (UpdateMeanFreePathForLatticeChangeover(aTrack)){
//    G4CMPVProcess::UpdateSCAfterLatticeChange();
//  }
//
  *condition = Forced;
  return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


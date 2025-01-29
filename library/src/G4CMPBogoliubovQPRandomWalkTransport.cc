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
  //fStepX =  0.0;
  //fStepY =  0.0;
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
  G4cout << "REL-- InRandomWalkTransport::StartTracking() A" << G4endl;
    G4VProcess::StartTracking(track);    // Apply base class actions
    LoadDataForTrack(track);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4CMPBogoliubovQPRandomWalkTransport::AlongStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double previousStepSize,
                             G4double currentMinimalStep,
                             G4double& currentSafety,
                             G4GPILSelection* selection)
{
  G4cout << "REL-- InRandomWalkTransport::AlongStepGPIL() A" << G4endl;

  G4cout << "---> REL/EY: At beginning of ASGPIL, track volume: " << track.GetVolume()->GetName() << G4endl;
  G4cout << "---> REL/EY: At beginning of ASGPIL, previousStepSize: " << previousStepSize << ", currentSafety: " << currentSafety << G4endl;
  G4cout << "---> REL/EY: At beginning of ASGPIL, momentum information: " << track.GetMomentumDirection() << G4endl;
  G4cout << "---> REL/EY: At beginning of ASGPIL, position: " << track.GetPosition() << G4endl;

  G4cout << "---> REL/EY: What is fGeomBoundary?" << fGeomBoundary << G4endl;
  
  //To determine if in turnaround step here, do the following:
  //1. Find current volume a la touchable.
  //2. Find momentum. If in turnaround step, current volume will be NOT the volume into which we're reflecting.
  //3. Find the volume that corresponds to a slight translation along the momentum direction. If it's NOT the same as the current
  //   volume, we're in a turnaround step and we set this to DBL_MAX for this step.
  //4. If else, we're not in a turnaround step -- we're potentially in a transmission step. Proceed as appropriate.
  //   -- Step is still not going to run okay here... we're passing in an INCOMPLETE step into the checkSurfaceNormal calculation. I think we need
  //      a CheckSurfaceNormal rendition for these stages where we need to calculate the 

  //Overall steps
  //0. Update the current lattice.
  //1. Determine if in turnaround step from reflection. If we are, set to DBL_MAX and skip
  //2. Otherwise, it's a transmission step.
  //3. Verify that the direction of momentum is going into a material with a lattice
  //4. Identify a stepX and step Y in that direction, making sure to draw Dn from the correct lattice

  //--------------------------------------------------------------------
  //0. Some preliminaries
  if(!fSafetyHelper) {
    G4TransportationManager* transportMgr ;
    transportMgr = G4TransportationManager::GetTransportationManager() ;
    fSafetyHelper = transportMgr->GetSafetyHelper();        
    fSafetyHelper->InitialiseHelper();
  }
  *selection = NotCandidateForSelection;
  fPathLength = currentMinimalStep;
  G4cout << "---> REL/EY: fPathLength in ASGPIL: " << currentMinimalStep << G4endl;

  // Do not need to check the material table since this is loaded in via
  // the logical and physical lattice
  // If lattice is a superconductor then fDn != 0 and fGapEnergy != 0
  G4double energy = track.GetKineticEnergy();
  
  //Get velocity of the current track
  G4double velocity = track.GetVelocity();
  
  //Initialize the lattice manager
  G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();

  //Update the lattice in procUtils here!

  
  //--------------------------------------------------------------------  
  //1. Determine if we're in turnaround step.
  G4VPhysicalVolume * currentVolume = track.GetVolume();
  G4ThreeVector trackPosition = track.GetPosition();
  G4ThreeVector momentumDir = track.GetMomentumDirection();
  G4ThreeVector trackPosition_eps = trackPosition+momentumDir*0.01*nm;
  G4ThreeVector trackPosition_mineps = trackPosition-momentumDir*0.01*nm; 
  G4VPhysicalVolume * currentVolPlusEps = G4CMP::GetVolumeAtPoint(trackPosition_eps);
  G4VPhysicalVolume * currentVolMinEps = G4CMP::GetVolumeAtPoint(trackPosition_mineps);
  G4StepStatus theStatus = track.GetStep()->GetPreStepPoint()->GetStepStatus();
  G4cout << "---> REL/EY: CurrentVolume by track: " << currentVolume->GetName() << ", CurrentVolumePlusEps, by GetVolumeAtPoint: " << currentVolPlusEps->GetName() << ", currentVolumeMinusEps, by GetVolumeAtPoint: " << currentVolMinEps->GetName() << ",  step status: " << track.GetStep()->GetPreStepPoint()->GetStepStatus() << G4endl;

  
  
  //If these are not the same, then we're in a turnaround step. This assumes that the distance to the next geometric feature are not smaller than 0.1 nm.
  fTrackOnBoundary = false;
  if( currentVolPlusEps != currentVolume && theStatus == fGeomBoundary ){    
    G4cout << "---> REL/EY: In a turnaround step. Killing the transport GPIL." << G4endl;
    fTrackOnBoundary = true;
    //fStepX = 0.0;
    //fStepY = 0.0;
    isActive = false;
    return DBL_MAX;
  }
  if( theStatus == fGeomBoundary ) fTrackOnBoundary = true;
  
  G4cout << "---> REL/EY: Not in a turnaround step. Continuing the transport GPIL." << G4endl;

  G4cout << "---> REL/EY: theStatus: " << theStatus << G4endl;
  //--------------------------------------------------------------------  
  //2. Verify that, if we're on a boundary, the direction of momentum is going into a volume with a lattice. Since we've updated the
  //   lattice at this point, we should confirm that the current procUtils lattice is not null, and that it is the
  //   same as the lattice corresponding to what the latticeManager sees as belonging to currentVolPlusEps.
  //   This is purely a check.
  if( theStatus == fGeomBoundary && !LM->HasLattice(currentVolPlusEps) ){ //NEED TO DO THE ADDITIONAL CHECK ONCE WE UPDATE THE LATTICE FIRST
    G4ExceptionDescription msg;
    msg << "We're on a boundary (i.e. our boundary is now behind us) and find no lattice in the region where a BogoliubovQP is intending to go during the QP random walk transport's AlongStepGPIL function.";
    G4Exception("G4CMPBogoliubovQPRandomWalkTransport::AlongStepGetPhysicalInteractionLength", "BogoliubovQPRandomWalkTransport001",FatalException, msg);
  }

 
  //3. Identify how best to access the diffusion constant. For now, let's access just from the lattice manager, and not from a procUtils
  //   stored lattice info. If we're on a geometry boundary, then do use the shifted volume. Otherwise, use the currentVolume a la the touchables.
  G4LatticePhysical* theLat;
  if( theStatus == fGeomBoundary ){
    //Access the lattice parameters using the shifted volume
    theLat = LM->GetLattice(currentVolPlusEps);
    fGapEnergy = theLat->GetSCDelta0() * 0.98; //REL TEMPORARY TESTING KLUDGE
    fTau0_qp = theLat->GetSCTau0qp();
    fDn = theLat->GetSCDn();
  }
  else{
    //Access the lattice parameters using the touchable volume
    theLat = LM->GetLattice(currentVolume);
    fGapEnergy = theLat->GetSCDelta0() * 0.98; //REL TEMPORARY TESTING KLUDGE
    fTau0_qp = theLat->GetSCTau0qp();
    fDn = theLat->GetSCDn();
  }

  G4cout << "---> REL/EY: Gap energy: " << fGapEnergy << ", energy: " << energy << ", DN: " << fDn << G4endl; 
  
  //4. Given this, compute a step length based on the time step.
  if ((energy>=fGapEnergy) && (fDn)>0){ 
    isActive = true;
    *selection = NotCandidateForSelection;
    
    //Calculate the energy dependent diffusion constant
    G4double E_ratio = fGapEnergy/energy;
    G4double E_ratio2 = pow(E_ratio,2.0);
    fDiffConst = fDn*sqrt(1-E_ratio2);
    G4double E_ratio_diff = energy/fGapEnergy - 1;
    fTimeStep = fPathLength/velocity;
    
    
    //Now we can calculate the std of the guassian for the RW.
    //Note that here, we don't need to restrict this yet. If this is going to hit a boundary, we'll adjust the time step (i.e. velocity) during
    //the AlongStepDoIt function.
    G4double sigma = sqrt(2.0*fDiffConst*fTimeStep);
    G4cout << "--->REL/EY: fTimeStep: " << fTimeStep << ", velocity: " << velocity << ", sigma: " << sigma << ", energy: " << energy << ", fDn: " << fDn << G4endl;
    G4RandGauss* gauss_dist = new G4RandGauss(G4Random::getTheEngine(),0.0,sigma);
    fPathLength = fabs(gauss_dist->fire());
    G4cout << "--->REL/EY: Successfully returning AlongStepGPIL path length of " << fPathLength << G4endl;
    G4cout << "--->REL/EY: Momentum direction is: " << momentumDir.unit() << G4endl;    
    return fPathLength;  
  }
  else{
    G4cout << "---> REL/EY: QP energy is too low or we're missing a Dn." << G4endl;
    //fStepX = 0.0;
    //fStepY = 0.0;
    isActive = false;
    return DBL_MAX;
  }
}

    
    /*
    fStepX = gauss_dist->fire();
    fStepY = gauss_dist->fire();
    fPathLength = sqrt((fStepX*fStepX+fStepY*fStepY));
    
    G4cout << "--->REL/EY: fStepX: " << fStepX << ", fStepY: " << fStepY << ", fPathLength: " << fPathLength << G4endl;
    
    //At this point, we're not in a turnaround step (otherwise we would have returned above). This means that if we're on a boundary,
    //we're either in a transmission step or in the post-turnaround step involved in reflection. However, in both of
    //those cases, we're *in* the volume into which we're moving. So at THIS point we need to "look back" to the
    //volume behind us to make sure we have both points needed to find the surface norm.    
    if( theStatus == fGeomBoundary ){

      G4cout << "---> REL/EY: Computing the surface norm for a guessed boundary behind us." << G4endl;
      
      //Calculate a surface norm using a single pre-step point and a "guess" at the direction of the surface behind us
      G4ThreeVector surfaceNorm = G4CMP::GetSurfaceNormal(track.GetStep()->GetPreStepPoint(),-1*momentumDir);
      G4cout << "---> REL/EY: Surface norm: " << surfaceNorm << G4endl;

      //Now we need to make sure our sampled point is in the same direction (relative to the surface) as the initial velocity
      //given to us by the boundary process. The way to think about this is using three vectors:
      //1. The surface normal
      //2. The momentum direction.
      //3. The vector of new sampled position.
      //Since the surface normal may point either into or out of a volume, and since we can nest volumes in confusing ways, the
      //important thing to recognize is the *relationship between the surface normal and the momentum direction*.

      //If the momentum direction and surface normal are in the same half-plane (I.e. have positive dot product), then we will need the
      //sampled position vector to also have positive dot product with the surface normal.
      double surfaceNormalMomentumDirDotProduct = surfaceNorm.dot(momentumDir);
      if( surfaceNormalMomentumDirDotProduct > 0 ){
	G4ThreeVector nextStepVect(fStepX,fStepY,0);
	do{
	  G4double phi = 2*pi*G4UniformRand();
	  fStepX = fPathLength*std::cos(phi);
	  fStepY = fPathLength*std::sin(phi);
	  nextStepVect.setX(fStepX);
	  nextStepVect.setY(fStepY);
	}
	while( surfaceNorm.dot(nextStepVect) <= 0 );
      }

      //If the momentum direction and surface normal are in opposite half-planes (i.e. have negative dot product), then we will need the
      //sampled position vector to also have negative dot product with the surface normal.
      else if( surfaceNormalMomentumDirDotProduct < 0 ){
	G4ThreeVector nextStepVect(fStepX,fStepY,0);
	do{
	  G4double phi = 2*pi*G4UniformRand();
	  fStepX = fPathLength*std::cos(phi);
	  fStepY = fPathLength*std::sin(phi);
	  nextStepVect.setX(fStepX);
	  nextStepVect.setY(fStepY);
	}
	while( surfaceNorm.dot(nextStepVect) >= 0 );
      }
      //They're parallel? What?
      else{
	G4ExceptionDescription msg;
	msg << "Somehow the surface normal dotted into the momentum direction is zero, which means something is wrong.";
	G4Exception("G4CMPBogoliubovQPRandomWalkTransport::AlongStepGetPhysicalInteractionLength", "BogoliubovQPRandomWalkTransport002",FatalException, msg);
      }
    }


    //If we're not on a boundary, then we can just let our X and Y steps be what they were computed to be earlier. We should also try to update
    //the direction information. Let's see if this works...
    double dX = fStepX;
    double dY = fStepY;
    G4ThreeVector newDir(dX,dY,0);    
    //    fParticleChange.ProposeMomentumDirection(newDir.unit());
    fSafetyHelper->Locate(trackPosition,newDir.unit()); //REL being a cowboy here... will this work?
    G4cout << "--->REL/EY: Attempting a fSafetyHelper->Locate() with track position: " << trackPosition << " and orientation: " << newDir.unit() << G4endl;
    */
/*
    
    G4bool positionPlusEpsLattice =LM->HasLattice(currentVol_eps);
    G4bool positionMinusEpsLattice =LM->HasLattice(currentVol_minus_eps);

    G4cout << "current volume is: " << currentVol->GetName() << " current volume eps is: " << currentVol_eps->GetName() << G4endl;
    //if(verboseLevel>2) {
        G4cout << "current position has lattice status "<<LM->HasLattice(currentVol)<<G4endl;
        G4cout << "current position plus eps has lattice status "<<positionPlusEpsLattice<<G4endl;
        G4cout << "current position minus eps has lattice status "<<positionMinusEpsLattice<<G4endl;
	//}
    
    //Check if the current step is on the boundary
    //it is likely that the below boolean can just check if the current steppoint has fGeomBoundary status
    fTrackOnBoundary = currentVol_eps != currentVol_minus_eps;
    G4cout << "current volume, 2 is: " << currentVol->GetName() << ", current volume eps is: " << currentVol_eps->GetName() << G4endl;
    
    G4ThreeVector surfaceNorm;

    G4cout << "---> REL RWT A" << G4endl;

    G4cout << "Current minimal step: " << currentMinimalStep << G4endl;
    G4cout << "Track step pre-touchable volume: " << track.GetStep()->GetPreStepPoint()->GetTouchable()->GetVolume()->GetName() << G4endl;
    G4cout << "Track step post-touchable volume: " << track.GetStep()->GetPostStepPoint()->GetTouchable()->GetVolume()->GetName() << G4endl;
    G4cout << "Track step prepoint 2: " << track.GetStep()->GetPreStepPoint()->GetPhysicalVolume()->GetName() << G4endl;
    G4cout << "Track step postpoint 2: " << track.GetStep()->GetPostStepPoint()->GetPhysicalVolume()->GetName() << G4endl;
    
    G4cout << "---> REL RWT B" << G4endl;
   
    G4double proposedTimeStep = 0.0;
    
    G4LatticePhysical* theLat;

    G4cout << "---> REL track.prestepPoint(): " << track.GetStep()->GetPreStepPoint()->GetPosition() << ", postStepPoint: " << track.GetStep()->GetPostStepPoint()->GetPosition() << G4endl;
    G4ThreeVector preSP = track.GetStep()->GetPreStepPoint()->GetPosition();
    G4ThreeVector postSP = track.GetStep()->GetPostStepPoint()->GetPosition();
    G4ThreeVector disp = preSP - postSP;
    double realStepLength = disp.mag();
    G4cout << "---> REL track.stepLength(): " << track.GetStep()->GetStepLength() << ", mag of subtraction: " << realStepLength << G4endl;

    if (LM->HasLattice(currentVol) ){ 
        //Current position volume has a lattice attached
        theLat = LM->GetLattice(currentVol);
        fGapEnergy = theLat->GetSCDelta0();
        fTau0_qp = theLat->GetSCTau0qp();
        fDn = theLat->GetSCDn();
        proposedTimeStep = 0.5*fTau0_qp;
    } else{
        //The current volume doesnt have a lattice attached return max double and make process inactive
        //if(verboseLevel>2)
	  G4cout << "G4CMPBogoliubovQPRandomWalkTransport: inside AlongStepGPIL " <<"currentVol does not have a lattice therefore making along step inactive " << G4endl;
        fStepX = 0.0;
        fStepY = 0.0;
        isActive = false;
        return DBL_MAX;
    }

    //REL was before previous block a bit ago
    if (fTrackOnBoundary){
        surfaceNorm=G4CMP::GetSurfaceNormal(*track.GetStep());
        //if(verboseLevel>2) {
            G4cout << "G4CMPBogoliubovQPRandomWalkTransport: inside AlongStepGPIL : fTrackOnBoundary is true " <<"surfaceNorm : "<<surfaceNorm << G4endl;
            G4cout << "G4CMPBogoliubovQPRandomWalkTransport: inside AlongStepGPIL : Step status fGeomBoundary? : " << (track.GetStep()->GetPreStepPoint()->GetStepStatus()==fGeomBoundary) << G4endl;
            G4cout << "G4CMPBogoliubovQPRandomWalkTransport: inside AlongStepGPIL " <<"momentum : "<<pDir << G4endl;
            G4cout << "G4CMPBogoliubovQPRandomWalkTransport: inside AlongStepGPIL " <<"momentum dot surfaceNorm : "<<pDir.dot(surfaceNorm) << G4endl;
	    //}
    }

    
    G4cout << "---> REL RWT C" << G4endl;
    
    if ((energy>=fGapEnergy) && (fDn)>0){
        isActive = true;
        //*selection = CandidateForSelection;
	*selection = NotCandidateForSelection;
        
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
        //fTimeStep = std::min(proposedTimeStep, fPathLength/velocity);
	fTimeStep = fPathLength/velocity;
        
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

    G4cout << "---> REL RWT D" << G4endl;


}
*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4CMPBogoliubovQPRandomWalkTransport::PostStepGetPhysicalInteractionLength(
              const G4Track&, G4double, G4ForceCondition* condition)
{
  G4cout << "REL-- InRandomWalkTransport::PostStepGetPhysicalInteractionLength() A" << G4endl;
  *condition = NotForced;
  return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4CMPBogoliubovQPRandomWalkTransport::AlongStepDoIt(const G4Track& track, const G4Step& step) {

  G4cout << "REL-- InRandomWalkTransport::AlongStepDoIt() A" << G4endl;
  G4cout << "REL: firstlook status1: " << track.GetStep()->GetPreStepPoint()->GetStepStatus() << G4endl;
  G4cout << "REL: firstlook status2: " << step.GetPreStepPoint()->GetStepStatus() << G4endl;
  G4cout << "REL: firstlook status3: " << step.GetPostStepPoint()->GetStepStatus() << G4endl;


  
  fParticleChange.ProposeMomentumDirection(
    step.GetPostStepPoint()->GetMomentumDirection());
  fNewPosition = step.GetPostStepPoint()->GetPosition();
  fNewDirection = step.GetPostStepPoint()->GetMomentumDirection();
  fParticleChange.ProposePosition(fNewPosition);
  G4double velocity = step.GetPostStepPoint()->GetVelocity();
  fParticleChange.ProposeVelocity(velocity);
  fPositionChanged = false;

  G4double stepLength = step.GetStepLength();
  G4double epsilon = 1.0*nm;
  
  // Check if the particle met conditions to do random walk from GPIL command
  if(!isActive) {
    fPathLength = stepLength;
    G4cout << "Step is not active" << G4endl;
    G4cout << "Particle path length: " << fPathLength <<G4endl;
    G4cout << "Particle velocity: " << velocity <<G4endl;
    G4cout << "Time Step : " << fPathLength/velocity <<G4endl;
    G4cout << "Particle direction: " << fNewDirection << G4endl;
  }

  // Particle did meet conditions to undergo RW
  else {

    //Update the velocity
    velocity = fPathLength / fTimeStep;
    fParticleChange.ProposeVelocity(velocity);

    //Need to check to see whether our proposed step is within the bounds of our current volume. fPathLength is passed as a const.
    fOldPosition = step.GetPreStepPoint()->GetPosition();
    fNewDirection = step.GetPreStepPoint()->GetMomentumDirection();
    fNewPosition = fOldPosition+fPathLength*fNewDirection;
    G4double dist = 0.0;
    G4double safety = fSafetyHelper->ComputeSafety(fOldPosition);

    //Test. I think CheckNextStep may also be valuable here, especially when we get into scenarios where we enter daughter volumes
    double nextStepSafety = 0;
    double nextStepLength = fSafetyHelper->CheckNextStep(fOldPosition,fNewDirection,fPathLength,nextStepSafety);
    G4cout << "---> REL/EY: Checking next step. NextStepLength: " << nextStepLength << ", nextStepSafety: " << nextStepSafety << G4endl;

    //We're actually going to try using checkNextStep instead. Here, if we're not on a boundary, the returned step should be kInfinity and
    //the safety should be nonzero.
    if( nextStepLength == kInfinity ){ //We're out in the bulk
      G4cout << "---> REL/EY: the proposed step does not cross a boundary. Setting position manually." << G4endl;
      fSafetyHelper->ReLocateWithinVolume(fNewPosition);
      fParticleChange.ProposeMomentumDirection(fNewDirection);
      fPositionChanged = true;
    }
    else if( nextStepLength != kInfinity ){ //We're on a surface
      if( step.GetPostStepPoint()->GetStepStatus() != fGeomBoundary ){
	G4ExceptionDescription msg;
	msg << "Somehow the CheckNextStep returned a step length that is not kInfinity but the step is indeed boundary limited. Are we actually landing on a boundary?";
	G4Exception("G4CMPBogoliubovQPRandomWalkTransport::AlongStepGetPhysicalInteractionLength", "BogoliubovQPRandomWalkTransport00X",FatalException, msg);
      }
      else{
	G4cout << "---> REL/EY: CheckNextStep shows that we've hit a boundary with our fPathLength computed in G4CMPBogoliubovRandomWalkTransport class. Recalculating time-to-wall and required velocity for this, but not setting new position -- will let Transportation set this." << G4endl;
	//Compute a representative deltaT at which we reach the wall given diffusion. Note that here, we really should be drawing from a
	//distribution, since there's not a one-to-one for the time at which diffusion from one point to another occurs. But for now we'll
	//make it a single time for simplicity.
	G4double newTime = dist*dist / 2.0 / fDiffConst;
	velocity = dist / newTime;
	fParticleChange.ProposeVelocity(velocity);
	fParticleChange.ProposeMomentumDirection(fNewDirection);
	//Here, we do not propose a new position -- we let transport do its thing, given that it has ostensibly won the GPIL race.
      }
    }
    
    /*
    if( fSafetyHelper->RecheckDistanceToCurrentBoundary(fOldPosition, fNewDirection, fPathLength, &dist, &safety) ){

      G4cout << "---> REL/EY: For an old position of: " << fOldPosition << " and a new direction of " << fNewDirection << ", the  ReCheckDistanceToCurrentBoundary yields a dist of " << dist << G4endl;
      
      //If the randomly drawn diffusion path length is less than the distance to the wall, then we just drop our pin and call it a day
      if (fPathLength < dist){
	G4cout << "---> REL/EY: fPathLength is LESS THAN the distance to wall along the momentum direction. Setting position manually." << G4endl;
	fSafetyHelper->ReLocateWithinVolume(fNewPosition);
	fParticleChange.ProposeMomentumDirection(fNewDirection);
	fPositionChanged = true;
      }
      //If the randomly drawn diffusion path length is equal to the distance to the wall, then we should allow ourselves to run whatever
      //post-step process is limiting this step. If a phonon needs to be created, let's let one be created. For now, let's just relocate to
      //ALMOST the wall, and make the fPathLength just a bit shorter so that it for sure beats out Transportation.
      else if (fPathLength == dist ){
	G4cout << "---> REL/EY: fPathLength is EQUAL to the distance to wall along the momentum direction." << G4endl;
	G4ThreeVector amendedPosition = fNewPosition - epsilon * fNewDirection;
	fNewPosition = amendedPosition;
	fSafetyHelper->ReLocateWithinVolume(fNewPosition);
	fParticleChange.ProposeMomentumDirection(fNewDirection);
	fPathLength -= epsilon;
	fPositionChanged = true;
      }
      //If the randomly drawn diffusion path length is larger than the distance to the wall, then we need to ensure that transportation wins.
      //Here, we don't modify the particlechange but we do need to update the velocity and the time step so that when Transportation wins, it
      //"gets to the boundary" at the appropriate time given diffusion.
      else{
	G4cout << "---> REL/EY: fPathLength is LARGER than the distance to the wall along the momentum direction. Recalculating time-to-wall and required velocity for this, but not setting new position -- will let Transportation set this." << G4endl;
	//Compute a representative deltaT at which we reach the wall given diffusion. Note that here, we really should be drawing from a
	//distribution, since there's not a one-to-one for the time at which diffusion from one point to another occurs. But for now we'll
	//make it a single time for simplicity.
	G4double newTime = dist*dist / 2.0 / fDiffConst;
	velocity = dist / newTime;
	fParticleChange.ProposeVelocity(velocity);
	fParticleChange.ProposeMomentumDirection(fNewDirection);
	//Here, we do not propose a new position -- we let transport do its thing, given that it has ostensibly won the GPIL race.
      }
    }
    //If we fail the safety check... what happens?
    else{
      //For now, just yell
      G4ExceptionDescription msg;
      msg << "We've somehow failed the step of checking if fPathLength is consistent with remaining in our current volume.";
      G4Exception("G4CMPBogoliubovQPRandomWalkTransport::AlongStepGetPhysicalInteractionLength", "BogoliubovQPRandomWalkTransport00X",FatalException, msg);

    }
    */
    
  }

  //If we've manually changed the position in such a way that transportation doesn't win the GPIL race, then propose the new position here.
  if(fPositionChanged){
    fParticleChange.ProposePosition(fNewPosition);
  }

  G4cout << "---> REL/EY: At end of AlongStepDoIt, old position: " << fOldPosition << G4endl;
  G4cout << "---> REL/EY: At end of AlongStepDoIt, new position: " << fNewPosition << G4endl;
  
  //Should this go here or in the above block?
  fParticleChange.ProposeTrueStepLength(fPathLength);
  G4cout << "---> REL/EY: Proposing true step length of fPathLength: " << fPathLength << " in the direction of " << fNewDirection << G4endl;
  
  return &fParticleChange;
}      

      
  /*
      

    
    // Update relevant particle changes
    //'Velocity' of particle is set by displacement and simulation time step
    G4cout << "--->REL/EY, fTimeStep: " << fTimeStep << G4endl;
    velocity = fPathLength/fTimeStep;
    fParticleChange.ProposeVelocity(velocity);
      
    //'MomentumDirection' of particle is set by the random walk gaussian steps
    fNewDirection = G4ThreeVector(fStepX/fPathLength,fStepY/fPathLength,0.0);
    fOldPosition = step.GetPreStepPoint()->GetPosition();
    fNewPosition = fOldPosition+fPathLength*fNewDirection;
      
    G4double dist = 0.0;
    G4double safety = fSafetyHelper->ComputeSafety(fOldPosition);

    G4cout << "REL: fOldPosition, prior to rechecks: " << fOldPosition << G4endl;
    G4cout << "REL: fNewPosition, prior to rechecks: " << fNewPosition << G4endl;
    G4cout << "REL: safety:  " << safety << G4endl;
    G4cout << "REL: status1: " << track.GetStep()->GetPreStepPoint()->GetStepStatus() << G4endl;
    G4cout << "REL: status2: " << step.GetPreStepPoint()->GetStepStatus() << G4endl;
    G4cout << "REL: status3: " << step.GetPostStepPoint()->GetStepStatus() << G4endl;

    
    //    if(fSafetyHelper->RecheckDistanceToCurrentBoundary(fOldPosition, fNewDirection, fPathLength, &dist, &safety)){
    if(fSafetyHelper->RecheckDistanceToCurrentBoundary(fOldPosition, fNewDirection, fPathLength, &dist, &safety)){
      if (fPathLength <=dist){
	G4cout << "REL fPathLength, " << fPathLength << ", is less than or equal to dist: " << dist << G4endl;
	G4cout << "REL new safety, after the rechecking: " << safety << G4endl;
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
	  G4cout << "REL Redoing safety check: " << fSafetyHelper->ComputeSafety(fNewPosition) << G4endl;
	}
      }else{
	//if(verboseLevel>2){
	G4cout << "G4CMPBogoliubovQPRandomWalkTransport: inside AlongStepDoIt, step colliding with boundary  \t dist : "<<dist << G4endl;
	G4cout << "\t new position from GPIL : "<< fNewPosition << G4endl;
	G4cout << "\t proposed new position : "<< fOldPosition+(dist)*fNewDirection<<G4endl;
	// }
	fNewPosition = fOldPosition+(dist)*fNewDirection;
	//fSafetyHelper->ReLocateWithinVolume(fNewPosition); //REL added -- should this be here?
	if(verboseLevel>2) G4cout << "\t new position : "<< fNewPosition << G4endl;
	fParticleChange.ProposeMomentumDirection(fNewDirection);
	fPositionChanged = true;
	fPathLength=dist;
      }
      
    }

    
    
    //if(verboseLevel>2){
    G4cout << "G4CMPBogoliubovQPRandomWalkTransport: inside AlongStepDoIt  \t Old position : "<< fOldPosition <<G4endl;
    G4cout << "\t New position : "<< fNewPosition <<G4endl;
    G4cout << "\t New safety : "<< safety <<G4endl;
       ///}
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
  */

  
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
  G4cout << "REL-- InRandomWalkTransport::GetContinuousStepLimit() A" << G4endl;
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


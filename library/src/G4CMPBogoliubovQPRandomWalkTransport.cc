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
  : G4VContinuousDiscreteProcess(name,fPhonon),G4CMPBoundaryUtils(this),G4CMPSCUtils(),
  fNewPosition(0.,0.,0.),
  fNewDirection(1.,0.,0.)
{
  verboseLevel = G4CMPConfigManager::GetVerboseLevel();
  SetProcessSubType(fBogoliubovQPRandomWalkTransport);
    
  //This time step will be overwritten by the step-limiting length (discrete process GPIL race)
  fTimeStep = 0;
  fPathLength =  0.0;
  fPreDiffusionPathLength = 0.0;
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
  G4cout << "---> REL/EY: At beginning of ASGPIL, global time: " << track.GetGlobalTime() << G4endl;
  
  
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

  //Set the path length and pre-diffusion path length to the length suggested by discrete process race winner
  fPathLength = currentMinimalStep;
  fPreDiffusionPathLength = currentMinimalStep;
  G4cout << "---> REL/EY: fPathLength in ASGPIL: " << currentMinimalStep << G4endl;
  G4cout << "---> REL/EY: fPreDiffusionPathLength in ASGPIL: " << currentMinimalStep << G4endl;
  G4cout << "---> REL/EY: velocity in ASGPIL: " << track.GetVelocity() << G4endl;

  
  //Get energy and velocity of track
  G4double energy = track.GetKineticEnergy();
  G4double velocity = track.GetVelocity();
  
  //Initialize the lattice manager
  G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();

  
  //--------------------------------------------------------------------  
  //1. Determine if we're in turnaround step.

  //REL 0.01 nm has some odd edge cases pop up where currentVolPlusEps isn't calculated properly. Guessing it's something in that function
  //that doesn't like the smaller tolerance because 1.0*nm works okay.
  G4double boundTolerance = 1.0*nm;//0.01*nm; 
  G4VPhysicalVolume * currentVolume = track.GetVolume();
  G4ThreeVector trackPosition = track.GetPosition();
  G4ThreeVector momentumDir = track.GetMomentumDirection();
  G4ThreeVector trackPosition_eps = trackPosition+momentumDir*boundTolerance;
  G4ThreeVector trackPosition_mineps = trackPosition-momentumDir*boundTolerance;
  G4VPhysicalVolume * currentVolPlusEps = G4CMP::GetVolumeAtPoint(trackPosition_eps);
  G4VPhysicalVolume * currentVolMinEps = G4CMP::GetVolumeAtPoint(trackPosition_mineps);
  G4StepStatus theStatus = track.GetStep()->GetPreStepPoint()->GetStepStatus();
  G4cout << "---> REL/EY: momentumDir in ASGPIL: " << momentumDir << G4endl;
  G4cout << "---> REL/EY: CurrentVolume by track: " << currentVolume->GetName() << ", CurrentVolumePlusEps, by GetVolumeAtPoint: " << currentVolPlusEps->GetName() << ", currentVolumeMinusEps, by GetVolumeAtPoint: " << currentVolMinEps->GetName() << ",  step status: " << track.GetStep()->GetPreStepPoint()->GetStepStatus() << G4endl;
  
  //If the volume plus epsilon (in the direction of the step) and the current volume are not the same (and if we're
  //in a boundary-limited step), then we're in a turnaround step. This assumes that the real distance to the next
  //geometric feature is not smaller than the boundTolerance -- in that case, this may break down a bit. Return
  //that this is inactive.
  fTrackOnBoundary = false;
  if( currentVolPlusEps != currentVolume && theStatus == fGeomBoundary ){    
    G4cout << "---> REL/EY: In a turnaround step. Killing the transport GPIL." << G4endl;
    fTrackOnBoundary = true;
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
  if( theStatus == fGeomBoundary && !LM->HasLattice(currentVolPlusEps) ){ 
    G4ExceptionDescription msg;
    msg << "We're on a boundary (i.e. our boundary is now behind us) and find no lattice in the region where a BogoliubovQP is intending to go during the QP random walk transport's AlongStepGPIL function.";
    G4Exception("G4CMPBogoliubovQPRandomWalkTransport::AlongStepGetPhysicalInteractionLength", "BogoliubovQPRandomWalkTransport001",FatalException, msg);
  }

  //3. The diffusion constant and other superconductor-specific information are accessed via the
  //   SCUtils class that this inherits from. No need to actually grab it from somewhere: it's already
  //   a data member.
  G4cout << "---> REL/EY: The gap energy, drawn from SCUtils, is: " << fGapEnergy << G4endl;
  G4cout << "---> REL/EY: The Dn, drawn from SCUtils, is: " << fDn << G4endl;
  G4cout << "---> REL/EY: The Teff, drawn from SCUtils, is: " << fTeff << G4endl;  
  G4cout << "---> REL/EY: Gap energy: " << fGapEnergy << ", energy: " << energy << ", DN: " << fDn << G4endl; 
  
  //4. Given this, compute a step length based on the time step.
  if ((energy>=fGapEnergy) && (fDn)>0){ 
    isActive = true;
    *selection = NotCandidateForSelection;
    
    //Calculate the energy dependent diffusion constant
    G4double E_ratio = fGapEnergy/energy;
    G4double E_ratio2 = pow(E_ratio,2.0);
    fDiffConst = fDn*sqrt(1-E_ratio2);
    fTimeStep = fPathLength/velocity;
    G4cout << "--->REL/EY: Diffusion constant (adjusted): " << fDiffConst << G4endl;

    //Now we can calculate the std of the gaussian for the RW. A few critical notes:
    //1. Since boundary interactions are handled between this function and the AlongStepDoIt (and critically rely on
    //   the Transportation process having its say), we don't restrict this distance here.
    //2. Since we're randomly drawing phi with our "pre-randomization" strategy, we're only sampling a 1D distribution
    //   in displacement. However, we can't just use a simple gaussian in r because it would over-weight low-r displacements.
    //   Here the density of points emerging in any direction from the pre-step points must fall off as a gaussian, which
    //   means that in radial coordinates our sampling function needs to pick up an additional factor of r. There are two
    //   ways to do this:
    //   a. Use a G4RandGeneral function and pass a Gauss * r function in (the "rigorous" way that is likely to be slower
    //      because we need a fine array for sampling)
    //   b. Sample in X and Y and use them to compute an R, ignoring the angle implied by X and Y. The resulting radial
    //      distribution should be the same as that in a, and this will give faster performance. This is what we will use.
    G4double sigma = sqrt(2.0*fDiffConst*fTimeStep);
    G4cout << "--->REL/EY: fTimeStep: " << fTimeStep << ", velocity: " << velocity << ", sigma: " << sigma << ", energy: " << energy << ", fDn: " << fDn << G4endl;
    G4RandGauss* gauss_dist = new G4RandGauss(G4Random::getTheEngine(),0.0,sigma);
    double gauss_dist_x = fabs(gauss_dist->fire());
    double gauss_dist_y = fabs(gauss_dist->fire());
    fPathLength = pow(gauss_dist_x*gauss_dist_x + gauss_dist_y*gauss_dist_y,0.5);
    G4cout << "--->REL/EY: Successfully returning AlongStepGPIL (diffusion-folded) path length of " << fPathLength << G4endl;
    G4cout << "--->REL/EY: Momentum direction is: " << momentumDir.unit() << G4endl;    
    return fPathLength;  
  }
  else{
    G4ExceptionDescription msg;
    msg << "QP energy is too low or we're missing a Dn. Returning DBL_MAX for alongStepDoIt.";
    G4Exception("G4CMPBogoliubovQPRandomWalkTransport::AlongStepGetPhysicalInteractionLength", "BogoliubovQPRandomWalkTransport002",FatalException, msg);
    isActive = false;
    return DBL_MAX;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4double 
G4CMPBogoliubovQPRandomWalkTransport::PostStepGetPhysicalInteractionLength(
              const G4Track& track, G4double previousStepSize, G4ForceCondition* condition)
{
  G4cout << "REL-- InRandomWalkTransport::PostStepGetPhysicalInteractionLength() A" << G4endl;

  //Since we're overriding this function as well, we'll have to call the MFP one (since it's called here in the base class).
  //This is all just to get the SCUtils updated at the beginning of the step calculus. Doing it in MFP because that's where
  //it's "usually" done with all of the other classes (even though it doesn't really matter which step it's called in here)
  double mfp = GetMeanFreePath(track, previousStepSize, condition);    
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

  //Get the initial and final times -- since Transport runs first in the alongStepDoIt processes,
  //we should have the "transport-limited" times
  G4double stepStartGlobalTime = step.GetPreStepPoint()->GetGlobalTime();
  G4double stepEndGlobalTime = step.GetPostStepPoint()->GetGlobalTime();
  G4double stepTransportOnlyDeltaT = stepEndGlobalTime-stepStartGlobalTime; //Limited by transport, needs to be overwritten.
  G4cout << "REL: stepStartGlobalTime: " << stepStartGlobalTime << ", stepEndGlobalTime: " << stepEndGlobalTime << ", stepTransportationOnlyDeltaT: " << stepTransportOnlyDeltaT << G4endl;

  G4ThreeVector preStepPoint = step.GetPreStepPoint()->GetPosition();
  G4ThreeVector postStepPoint = step.GetPostStepPoint()->GetPosition();
  G4double stepTransportationOnlyDeltaR = (postStepPoint-preStepPoint).mag();
  
  
  //This is the "old" velocity (i.e. that not computed using the diffusion step)
  G4double velocity = step.GetPostStepPoint()->GetVelocity();
  fParticleChange.ProposeVelocity(velocity);
  fPositionChanged = false;

  G4double stepLength = step.GetStepLength();
  G4double epsilon = 1.0*nm;
  
  // Check if the particle met conditions to do random walk from GPIL command. Here, this occurs
  // during a turnaround step where we want to set step lengths to zero. If we don't set the proposedTrueStepLength
  // to zero, then we'll end up mucking up our boundary processes which expect a zero step length
  // during turnaround steps.
  if(!isActive) {
    fPathLength = stepLength;
    fPreDiffusionPathLength = stepLength; 
    G4cout << "Step is not active" << G4endl;
    G4cout << "Particle path length: " << fPathLength <<G4endl;
    G4cout << "Particle pre-diffusion path length: " << fPreDiffusionPathLength <<G4endl;
    G4cout << "Particle velocity: " << velocity <<G4endl;
    G4cout << "Time Step : " << fPathLength/velocity <<G4endl;
    G4cout << "Particle direction: " << fNewDirection << G4endl;
    fParticleChange.ProposeLocalTime(0);
    fParticleChange.ProposeTrueStepLength(fPreDiffusionPathLength);
  }

  // Particle did meet conditions to undergo RW
  else {

    //Update the velocity
    //velocity = fPathLength / fTimeStep;
    G4cout << "---> REL/EY: in AlongStepDoit, velocity is set and proposed to " << velocity << G4endl;
    fParticleChange.ProposeVelocity(velocity);

    //Need to check to see whether our proposed step is within the bounds of our current volume. fPathLength is passed as a const, and here
    //is the diffusion displacement.
    fOldPosition = step.GetPreStepPoint()->GetPosition();
    fNewDirection = step.GetPreStepPoint()->GetMomentumDirection();
    fNewPosition = fOldPosition+fPathLength*fNewDirection;
    //G4double safety = fSafetyHelper->ComputeSafety(fOldPosition);

    G4cout << "---> REL/EY: In AlongStepDoIt, time of pre-step point is: " << step.GetPreStepPoint()->GetGlobalTime() << G4endl;
    G4cout << "---> REL/EY: In AlongStepDoIt, time of post-step point is: " << step.GetPostStepPoint()->GetGlobalTime() << G4endl;
    
    //Test. I think CheckNextStep may also be valuable here, especially when we get into scenarios where we enter daughter volumes
    double nextStepSafety = 0;    
    double nextStepLength = fSafetyHelper->CheckNextStep(fOldPosition,fNewDirection,fPathLength,nextStepSafety);
    G4cout << "---> REL/EY: Checking next step. CheckNextStep's next step length (dist to boundary): " << nextStepLength << ", nextStepSafety: " << nextStepSafety << G4endl;

    //We're actually going to try using checkNextStep instead. Here, if we're not on a boundary, the returned step should be kInfinity and
    //the safety should be nonzero.
    if( nextStepLength == kInfinity ){ //We're out in the bulk
      G4cout << "---> REL/EY: the proposed step does not cross a boundary. Setting position manually." << G4endl;
      fSafetyHelper->ReLocateWithinVolume(fNewPosition);
      fParticleChange.ProposeMomentumDirection(fNewDirection);

      //Since we are forced to use the pre-step point's velocity for propagation of this step (and the pre-step point's velocity is
      //the one we started with), transportation will add a time corresponding to traveling the calculated path length at that velocity.
      //We should subtract that time off our final proposed time.
      G4cout << "---> REL/EY: fTimeStep: " << fTimeStep << ", timeChangefromTransportationOnly: " << stepTransportOnlyDeltaT << G4endl;
      fParticleChange.ProposeLocalTime(fTimeStep-stepTransportOnlyDeltaT);


      //I think this needs to be set to the *old, pre-diffusion* path length, since other processes that aren't the step-limiting one will
      //need to subtract off a distance. That distance basically needs to be velocity * deltaT. The current process that limits the step
      //(here, not transportation) will zero out its number of interaction lengths.
      G4cout << "---> REL/EY: Proposing a true (diffusion-UNfolded) step length of: " << fPreDiffusionPathLength << G4endl;
      fParticleChange.ProposeTrueStepLength(fPreDiffusionPathLength);
      fPositionChanged = true;
    }
    else if( nextStepLength != kInfinity ){ //We're on a surface
      if( step.GetPostStepPoint()->GetStepStatus() != fGeomBoundary ){
	G4ExceptionDescription msg;
	msg << "Somehow the CheckNextStep returned a step length that is not kInfinity but the step is indeed boundary limited. Are we actually landing on a boundary?";
	G4Exception("G4CMPBogoliubovQPRandomWalkTransport::AlongStepDoIt", "BogoliubovQPRandomWalkTransport00X",FatalException, msg);
      }
      else{
	G4cout << "---> REL/EY: CheckNextStep shows that we've hit a boundary with our fPathLength computed in G4CMPBogoliubovRandomWalkTransport class. Recalculating time-to-wall and required velocity for this, but not setting new position -- will let Transportation set this." << G4endl;
	//Compute a representative deltaT at which we reach the wall given diffusion. Note that here, we really should be drawing from a
	//distribution, since there's not a one-to-one for the time at which diffusion from one point to another occurs. But for now we'll
	//make it a single time for simplicity.
	G4double newDiffusionTime = nextStepLength*nextStepLength / 2.0 / fDiffConst;

	//Sanity check: this step length and the step deltaT should be consistent with the velocity
	G4cout << "---> REL/EY: newDiffusionTime: " << newDiffusionTime << ", stepTransportationOnlyDeltaR: " << stepTransportationOnlyDeltaR << G4endl;
	fParticleChange.ProposeLocalTime(newDiffusionTime-stepTransportOnlyDeltaT);
	fParticleChange.ProposeMomentumDirection(fNewDirection);

	//We now modify the true step length (which becomes the next stepLength in the step) to be whatever the new diffusion time is, divided by the
	//velocity. This is working backwards from our above action of finding the step time from the velocity and the "true" step length, i.e. the one
	//that's given by the physics of the post-step processes. This "abridged true step length" will remove the necessary fPathLength from the scatter
	//despite Transportations' naive idea of what the next step length is (i.e. an "as the crow flies" length to the boundary).
	G4double abridgedTrueStepLength = velocity * newDiffusionTime;
	fParticleChange.ProposeTrueStepLength(abridgedTrueStepLength);
	G4cout << "---> REL/EY: proposed abridged diffusion-UNfolded step length is " << abridgedTrueStepLength << ", which is less than the original diffusion-UNfolded step length, " << fPreDiffusionPathLength << G4endl;
      }
    }
  }

  //If we've manually changed the position in such a way that transportation doesn't win the GPIL race, then propose the new position here.
  if(fPositionChanged){
    fParticleChange.ProposePosition(fNewPosition);
  }

  G4cout << "---> REL/EY: At end of AlongStepDoIt, old position: " << fOldPosition << G4endl;
  G4cout << "---> REL/EY: At end of AlongStepDoIt, new position: " << fNewPosition << G4endl;
  
  //Should this go here or in the above block?  
  return &fParticleChange;
}      

      

  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* 
G4CMPBogoliubovQPRandomWalkTransport::PostStepDoIt(const G4Track& track, const G4Step&)
{
  fParticleChange.Initialize(track);
  ClearNumberOfInteractionLengthLeft();		// All processes should do this! REL added but not sure it belongs here?
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
              const G4Track& track, G4double, G4ForceCondition* condition)
{
  G4cout << "REL Here in GetMeanFreePath for the G4CMPBogoliubovQPRandomWalkTransport." << G4endl;
  
  //This needs to be done so that we can update the SCUtils information, and since GetMeanFreePath for the discrete bit of this process should (?)
  //run first, it will hopefully happen before the AlongStepGPIL runs and requests the gap.
  if (UpdateMeanFreePathForLatticeChangeover(track)){
    UpdateSCAfterLatticeChange();
  }
  
  *condition = Forced;
  return DBL_MAX;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//This is an amended version of the logic block in G4CMPVProcess that allows us to update the lattice in procUtils and in SCUtils when it changes.
//Since we don't use a rateModel or other things explicitly here, it's not the exact same code, which is why I'm okayish with havine a separate
//funciton in this class.
G4bool G4CMPBogoliubovQPRandomWalkTransport::UpdateMeanFreePathForLatticeChangeover(const G4Track& aTrack)
{
  G4cout << "REL HereA_G4CMPBogoliubovQPRandomWalkTransport: loading data for track after lattice changeover, process: " << this->GetProcessName() << G4endl;
  G4cout << "Here, track length: " << aTrack.GetTrackLength() << G4endl;
  G4cout << "Current lattice a la lattice manager: " << G4LatticeManager::GetLatticeManager()->GetLattice(aTrack.GetVolume()) << G4endl;
    
  //Always do a check to see if the current lattice stored in this process is equal to the one that represents
  //the volume that we're in. Note that we can't do this with the "GetLattice()" and "GetNextLattice()" calls
  //here because at this point in the step, the pre- and post-step points both point to the same volume. Since
  //GetMeanFreePath is run at the beginning, I think the point at which a boundary interaction is assessed comes
  //later (hence why we can use that info in PostStepDoIts but not here.) Adding a statement about track length here,
  //since it seems that when a particle spawns it doesn't necessarily trigger this block, and I think we want it to.
  if( (((this->theLattice) && G4LatticeManager::GetLatticeManager()->GetLattice(aTrack.GetVolume())) &&
       (this->theLattice != G4LatticeManager::GetLatticeManager()->GetLattice(aTrack.GetVolume()))) ||
      aTrack.GetTrackLength() == 0.0 ){
    
    G4cout << "--------> REL the step length associated with this is " << aTrack.GetStep()->GetStepLength() << G4endl;
    
    
    
    //REL noting that if physical lattices are not 1:1 with volumes, something may get broken here... Should check a scenario of segmented SC...    
    this->LoadDataForTrack(&aTrack);
    G4cout << "REL G4CMPBogoliubovQPRandomWalkTransport: Successfully changed over to a new lattice for process " << this->GetProcessName() << G4endl;
    return true;
    
  }
  G4cout << "REL G4CMPBogoliubovQPRandomWalkTransport: Did not successfully change over to a new lattice for process " << this->GetProcessName() << G4endl;
  return false;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//This is meant to update superconductor info for the process if we move into a new lattice. Custom to this process, which does not have an associated
//rate model.
void G4CMPBogoliubovQPRandomWalkTransport::UpdateSCAfterLatticeChange()
{
  G4cout << "REL HereC_G4CMPVProcess: updating SC after lattice change" << G4endl;
  
  //First, determine if the new lattice is a SC. If not, then set the SCUtils info to null for this process  
  if( (this->theLattice)->GetSCDelta0() <= 0 ){
    this->SetCurrentSCInfoToNull();
    return;
  }

  //If it is a SC, then we should update the SC information for the SC utils class within the base of this.
  //Also, handle the checking/updating of the lookup tables to be used for each SC.
  this->LoadLatticeInfoIntoSCUtils(this->theLattice);
}


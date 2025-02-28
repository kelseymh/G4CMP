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
#include "G4RandomDirection.hh"

//Temporary REL
#include <fstream>

#include <iostream>
#include <cmath>
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4CMPBogoliubovQPRandomWalkTransport::G4CMPBogoliubovQPRandomWalkTransport(const G4String& name , G4CMPProcessSubType fBogoliubovQPRandomWalkTransport)
  : G4VContinuousDiscreteProcess(name,fPhonon),G4CMPBoundaryUtils(this),G4CMPSCUtils(),
  fNewPosition(0.,0.,0.),
  fNewDirection(1.,0.,0.)
{
  verboseLevel = 3;//G4CMPConfigManager::GetVerboseLevel();
  SetProcessSubType(fBogoliubovQPRandomWalkTransport);
    
  //This time step will be overwritten by the step-limiting length (discrete process GPIL race)
  fTimeStep = 0;
  fPathLength =  0.0;
  fPreDiffusionPathLength = 0.0;
  fDiffConst =  0.0;
  fBoundaryFudgeFactor = 1.0001;
  
  
  //fSafetyHelper is initialized in AlongStepGPIL
  fSafetyHelper=nullptr;
  
  //Temporary REL
  //fOutfile.open("/Users/ryanlinehan/QSC/Sims/Geant4/scRebuild-build/RandomWalkSampledDimlessTimes.txt",std::ios::trunc);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4CMPBogoliubovQPRandomWalkTransport::~G4CMPBogoliubovQPRandomWalkTransport()
{
  if(verboseLevel>2) {
    G4cout << "G4CMPBogoliubovQPRandomWalkTransport destruct " << GetProcessName()
          << G4endl;
  }

  //fOutfile.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CMPBogoliubovQPRandomWalkTransport::StartTracking(G4Track* track)
{
  //Debugging
  if( verboseLevel > 2 ){
    G4cout << "---------- G4CMPBogoliubovQPRandomWalkTransport::StartTracking ----------" << G4endl;
    G4cout << "ST Function Point A | StartTracking." << G4endl;
  }
  
  G4VProcess::StartTracking(track);    // Apply base class actions
  LoadDataForTrack(track);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//This is the "stretch" version
G4double G4CMPBogoliubovQPRandomWalkTransport::AlongStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double previousStepSize,
                             G4double currentMinimalStep,
                             G4double& currentSafety,
                             G4GPILSelection* selection)
{

  //Debugging
  if(verboseLevel>2){
    G4cout << "---------- G4CMPBogoliubovQPRandomWalkTransport::AlongStepGetPhysicalInteractionLength ----------" << G4endl;
    G4cout << "ASGPIL Function Point A | track volume: " << track.GetVolume()->GetName() << G4endl;
    G4cout << "ASGPIL Function Point A | previousStepSize: " << previousStepSize << ", currentSafety: " << currentSafety << G4endl;
    G4cout << "ASGPIL Function Point A | momentum direction: " << track.GetMomentumDirection() << G4endl;
    G4cout << "ASGPIL Function Point A | position: " << track.GetPosition() << G4endl;
    G4cout << "ASGPIL Function Point A | global time: " << track.GetGlobalTime() << G4endl;
  }
    
  //At this point, the currentMinimalStep is the one that has won the discrete GPIL race. 
  //Set the path length and pre-diffusion path length to the length suggested by discrete process race winner
  fPathLength = currentMinimalStep;
  fPreDiffusionPathLength = currentMinimalStep;
  if( verboseLevel>2){
    G4cout << "ASGPIL Function Point B | fPathLength = fPreDiffusionPathLength = currentMinimalStep: " << fPathLength << G4endl;
    G4cout << "ASGPIL Function Point B | velocity: " << track.GetVelocity() << G4endl;
  }

  //If we're in a turnaround step, then kill
  if( isActive == false ){

    //Some debugging
    if( verboseLevel>2 ){
      G4cout << "ASGPIL Function Point C | In a turnaround step. Killing the transport GPIL." << G4endl;
    }
    return DBL_MAX;
  }
  
  //Get some track and SC properties. The diffusion constant and other superconductor-specific information are accessed via the
  //SCUtils class that this inherits from. No need to actually grab it from somewhere: it's already a data member
  G4double energy = track.GetKineticEnergy();
  G4double velocity = track.GetVelocity();
  G4ThreeVector momentumDir = track.GetMomentumDirection();  
  if( verboseLevel > 2 ){
    G4cout << "ASGPIL Function Point D | Gap energy (drawn from SCUtils): " << fGapEnergy << ", particle energy: " << energy << G4endl;
    G4cout << "ASGPIL Function Point D | Dn (drawn from SCUtils): " << fDn << ", and fDiffConst: " << fDiffConst << G4endl;
    G4cout << "ASGPIL Function Point D | Teff (drawn from SCUtils): " << fTeff << G4endl;
  }
    
  //If our energy is appropriate and we can see diffusion info, trigger the "meat" of this function
  if ((energy>=fGapEnergy) && (fDn)>0){ 
    isActive = true;
    *selection = NotCandidateForSelection;

    //Debugging
    if( verboseLevel > 2 ){      
      G4cout << "ASGPIL Function Point E | isActive is true. currentMinimalStep/velocity: " << currentMinimalStep/velocity << ", fTimeStepToBoundary: " << fTimeStepToBoundary << G4endl;
    }

    //Here, we need to return a diffusion-displacement ("diffusion-folded") path length. If the currentMinimalStep matches
    //the fTimeStepToBoundary, it means that the winning discrete GPIL process is THIS process's discrete bit. In this case,
    //the diffusion displacement distance is simple: it's just the distance to the boundary (up to the fudge factor that we use
    //to actually trigger a boundary interaction in G4), and the time is just the time step to the boundary.
    double timeTolerance = 1E-10; //For floating point errors    
    if( (fabs(currentMinimalStep/velocity - fTimeStepToBoundary) < timeTolerance) || fVerySmallStep ){
      fPathLength = f2DSafety*fBoundaryFudgeFactor; //Fudge factor so that we can trigger transportation when we are pointed at the boundary exactly
      fTimeStep = fTimeStepToBoundary;

      //Debugging
      if( verboseLevel > 2 ){
	G4cout << "ASGPIL Function Point F | Looks like the boundary-limited case applies here. Returning fPathLength = f2DSafety = " << fPathLength << G4endl;
      }
    }
    //If another process wins the discrete GPIL race, the above block won't trigger. Here we need to compute a radius corresponding
    //to that (smaller) elapsed time. However, the radius *must* be less than the boundary-limited radius because otherwise the
    //particle would have passed through that radius, implying that the boundary-limited scenario *should* have won. To make sure
    //this is therefore self-consistent, we will draw a radius from a gaussian corresponding to the (smaller) elapsed time and
    //require the drawn radius to be smaller than the boundary distance.
    //NB: this makes qualitative sense but may be a place where we bias things on accident. May want to think through this logic a bit more.
    else{      
      fTimeStep = currentMinimalStep / velocity;
      G4double sigma1D = pow(2*fDiffConst*fTimeStep,0.5);
      G4RandGauss* gauss_dist = new G4RandGauss(G4Random::getTheEngine(),0.0,sigma1D);
      do{
	double gauss_dist_x = fabs(gauss_dist->fire());
	double gauss_dist_y = fabs(gauss_dist->fire());
	fPathLength = pow(gauss_dist_x*gauss_dist_x + gauss_dist_y*gauss_dist_y,0.5);
      }
      while( fPathLength > f2DSafety );

      //Debugging
      if( verboseLevel > 2 ){
	G4cout << "ASGPIL Function Point G | Sigma1D: " << sigma1D << " from fDiffConst: " << fDiffConst << G4endl;
	G4cout << "ASGPIL Function Point G | Looks like a different discrete process wins the GPIL race. Returning fPathLength = " << fPathLength << G4endl;
      }
    }

    //Debugging
    if( verboseLevel > 2 ){
      G4cout << "ASGPIL Function Point H | Successfully returning AlongStepGPIL (diffusion-folded) fPathLength of " << fPathLength << G4endl;
      G4cout << "ASGPIL Function Point H | Momentum direction is: " << momentumDir.unit() << G4endl;
    }
    return fPathLength;
  }

  //If we don't have a reasonable energy given the gap, or if we're missing Dn.
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
  //Debugging
  if( verboseLevel > 2 ){
    G4cout << "---------- G4CMPBogoliubovQPRandomWalkTransport::PostStepGetPhysicalInteractionLength ----------" << G4endl;
    G4cout << "PSGPIL Function Point A | In PostStepGetPhysicalInteractionLength" << G4endl;
  }

  //Since we're overriding this function as well, we'll have to call the MFP one (since it's called here in the base class).
  //This is all just to get the SCUtils updated at the beginning of the step calculus. Doing it in MFP because that's where
  //it's "usually" done with all of the other classes (even though it doesn't really matter which step it's called in here)
  double mfp = GetMeanFreePath(track, previousStepSize, condition);
  return mfp;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// "Stretch" version
G4VParticleChange* G4CMPBogoliubovQPRandomWalkTransport::AlongStepDoIt(const G4Track& track, const G4Step& step) {

  //Debugging
  if( verboseLevel > 2 ){
    G4cout << "---------- G4CMPBogoliubovQPRandomWalkTransport::AlongStepDoIt ----------" << G4endl;
    G4cout << "ASDI Function Point A | Step Status from track pre-step point: " << track.GetStep()->GetPreStepPoint()->GetStepStatus() << G4endl;
    G4cout << "ASDI Function Point A | Step Status from step pre-step point: " << step.GetPreStepPoint()->GetStepStatus() << G4endl;
    G4cout << "ASDI Function Point A | Step Status from step post-step point: " << step.GetPostStepPoint()->GetStepStatus() << G4endl;
  }

  //Particle change initialization and setting
  fParticleChange.Initialize(track);  
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
  G4ThreeVector preStepPoint = step.GetPreStepPoint()->GetPosition();
  G4ThreeVector postStepPoint = step.GetPostStepPoint()->GetPosition();
  G4double stepTransportationOnlyDeltaR = (postStepPoint-preStepPoint).mag();
  
  //Debugging
  if( verboseLevel > 2 ){
    G4cout << "ASDI Function Point B | stepStartGlobalTime: " << stepStartGlobalTime << G4endl;
    G4cout << "ASDI Function Point B | stepEndGlobalTime  : " << stepEndGlobalTime << G4endl;
    G4cout << "ASDI Function Point B | stepTransportOnlyDeltaT: " << stepTransportOnlyDeltaT << G4endl;
  }
  
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

    //Debugging
    if( verboseLevel > 2 ){
      G4cout << "ASDI Function Point C | Step is not active" << G4endl;
      G4cout << "ASDI Function Point C | fPathLength: " << fPathLength << G4endl;
      G4cout << "ASDI Function Point C | diffusion-unfolded (i.e. pre-diffusion) path length: " << fPreDiffusionPathLength << G4endl;
      G4cout << "ASDI Function Point C | velocity: " << velocity << G4endl;
      G4cout << "ASDI Function Point C | fPathLength/Velocity: " << fPathLength/velocity << G4endl;
      G4cout << "ASDI Function Point C | fNewDirection: " << fNewDirection << G4endl;
    }
    fParticleChange.ProposeLocalTime(0);
    fParticleChange.ProposeTrueStepLength(fPreDiffusionPathLength);
  }

  // Particle did meet conditions to undergo RW
  else {

    //Debugging
    if( verboseLevel > 2 ){
      G4cout << "ASDI Function Point D | velocity is set and proposed to " << velocity << ", which should never change." << G4endl;
    }
    fParticleChange.ProposeVelocity(velocity);

    //Take the old position and new direction and make a new position using the fPathLength, which here is the diffusion-folded
    //"displacement" distance
    fOldPosition = step.GetPreStepPoint()->GetPosition();
    fNewDirection = step.GetPreStepPoint()->GetMomentumDirection();
    fNewPosition = fOldPosition+fPathLength*fNewDirection;

    //Debugging
    if( verboseLevel > 2 ){
      G4cout << "ASDI Function Point E | time of pre-step point is: " << step.GetPreStepPoint()->GetGlobalTime() << G4endl;
      G4cout << "ASDI Function Point E | time of post-step point is: " << step.GetPostStepPoint()->GetGlobalTime() << G4endl;
    }

    //Test. I think CheckNextStep may also be valuable here, especially when we get into scenarios where we enter daughter volumes
    double nextStepSafety = 0;    
    double nextStepLength = fSafetyHelper->CheckNextStep(fOldPosition,fNewDirection,fPathLength,nextStepSafety);

    //Debugging
    if( verboseLevel > 2 ){
      G4cout << "ASDI Function Point F | Checking next step. CheckNextStep's next step length: " << nextStepLength << ", nextStepSafety: " << nextStepSafety << G4endl;
    }

    //We're actually going to try using checkNextStep instead. Here, if we're not on a boundary, the returned step should be kInfinity and
    //the safety should be nonzero.
    if( nextStepLength == kInfinity ){ //We're out in the bulk

      //Debugging
      if( verboseLevel > 2 ){
	G4cout << "ASDI Function Point G | the proposed step does not cross a boundary. Setting position manually." << G4endl;
      }
      fSafetyHelper->ReLocateWithinVolume(fNewPosition);
      fParticleChange.ProposeMomentumDirection(fNewDirection);

      //Since we are forced to use the pre-step point's velocity for propagation of this step (and the pre-step point's velocity is
      //the one we started with), transportation will add a time corresponding to traveling the calculated path length at that velocity.
      //We should subtract that time off our final proposed time. First, debugging.
      if( verboseLevel > 2 ){
	G4cout << "ASDI Function Point H | fTimeStep: " << fTimeStep << ", timeChangefromTransportationOnly: " << stepTransportOnlyDeltaT << G4endl;
	G4cout << "ASDI Function Point H | timeChangeFromTransportationOnly: " << stepTransportOnlyDeltaT << G4endl;
      }
      fParticleChange.ProposeLocalTime(fTimeStep-stepTransportOnlyDeltaT);


      //I think this needs to be set to the step-limiting *old, pre-diffusion* path length, since other processes that aren't the step-limiting one will
      //need to subtract off a distance. That distance basically needs to be velocity * deltaT. The current process that limits the step
      //(here, not transportation) will zero out its number of interaction lengths. First, debugging
      if( verboseLevel > 2 ){
	G4cout << "ASDI Function Point I | Proposing a true (diffusion-UNfolded) step length of: " << fPreDiffusionPathLength << G4endl;
      }
      fParticleChange.ProposeTrueStepLength(fPreDiffusionPathLength);
      fPositionChanged = true;
    }
    
    //If we hit the wall, we should have set things up so that:
    //1. The discrete process leading to a wall hit is the boundary-limited one set here
    //2. When we hit the wall, we set the time of impact to be consistent with diffusion,
    //   so that the above code (without the "relocate" bit) still applies.
    else if( nextStepLength != kInfinity ){ //We're on a surface
      if( step.GetPostStepPoint()->GetStepStatus() != fGeomBoundary ){
	G4ExceptionDescription msg;
	msg << "Somehow the CheckNextStep returned a step length that is not kInfinity but the step is indeed boundary limited. Are we actually landing on a boundary?";
	G4Exception("G4CMPBogoliubovQPRandomWalkTransport::AlongStepDoIt", "BogoliubovQPRandomWalkTransport00X",FatalException, msg);
      }
      else{

	//Debugging
	if( verboseLevel > 2 ){
	  G4cout << "ASDI Function Point J | CheckNextStep shows that we've hit a boundary. Using timestep computed in GetMFP, letting Transportation do move of particle." << G4endl;
	  G4cout << "ASDI Function Point J | fTimeStep: " << fTimeStep << G4endl;
	  G4cout << "ASDI Function Point J | timeChangeFromTransportationOnly: " << stepTransportOnlyDeltaT << G4endl;
	  G4cout << "ASDI Function Point J | Proposing true (diffusion-UNfolded) step length of: " << fPreDiffusionPathLength << G4endl;
	}
	  
	//Since we are forced to use the pre-step point's velocity for propagation of this step (and the pre-step point's velocity is
	//the one we started with), transportation will add a time corresponding to traveling the calculated path length at that velocity.
	//We should subtract that time off our final proposed time.
	fParticleChange.ProposeLocalTime(fTimeStep-stepTransportOnlyDeltaT);

	//I think this needs to be set to the *old, pre-diffusion* path length, since other processes that aren't the step-limiting one will
	//need to subtract off a distance. That distance basically needs to be velocity * deltaT. The current process that limits the step
	//(here, not transportation) will zero out its number of interaction lengths.
	fParticleChange.ProposeTrueStepLength(fPreDiffusionPathLength);
	fParticleChange.ProposeMomentumDirection(fNewDirection);      
      }
    }
  }

  //If we've manually changed the position in such a way that transportation doesn't win the GPIL race, then propose the new position here.
  if(fPositionChanged){
    fParticleChange.ProposePosition(fNewPosition);
  }

  //Debugging
  if( verboseLevel > 2 ){
    G4cout << "ASDI Function Point K | At end of ASDI, old position: " << fOldPosition << G4endl;
    G4cout << "ASDI Function Point K | At end of ASDI, new position: " << fNewPosition << G4endl;
  }
  
  //Should this go here or in the above block?  
  return &fParticleChange;
}      

      

  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* 
G4CMPBogoliubovQPRandomWalkTransport::PostStepDoIt(const G4Track& track, const G4Step&)
{
  //Debugging
  if( verboseLevel > 2 ){
    G4cout << "---------- G4CMPBogoliubovQPRandomWalkTransport::PostStepDoIt ----------" << G4endl;
  }
  
  G4double epsilon = 1*CLHEP::um; //REL HARDCODED, FIX

  //Determine if we're on a boundary. A few scenarios:
  //1. This shouldn't run if we are landing on a boundary in the step (where Transportation will be the thing
  //   that runs PostStepDoIt, not this function
  //   gets set to true in the AlongStepGPIL
  //2. This shouldn't run if we are in a turnaround step
  //3. This *can* run if we are limited by this process while *starting* from a boundary, so 
  //Determine if we're on a boundary. I think (?) that this should never get run if we *land* on a boundary,
  //and during the turnaround step, I think this doesn't run anyway (?). It's either the step that landed on the boundary
  //(in which transportation wins the post-step race), or the turnaround step (where we shouldn't run PostStepDoIt anyway)
  if( fTrackOnBoundary && !isActive ){
    G4ExceptionDescription msg;
    msg << "In a turnaround step and postStepDoIt is running";
    G4Exception("G4CMPBogoliubovQPRandomWalkTransport::PostStepDoIt", "BogoliubovQPRandomWalkTransport00X",FatalException, msg);
  }
  
  //Determine if we end up close to a boundary using Get2DSafety
  G4double the2DSafety = G4CMP::Get2DSafety(track.GetStep()->GetPreStepPoint()->GetTouchable(),
					    track.GetPosition(),
					    track.GetMomentumDirection(),
					    false);
  //Debugging
  if( verboseLevel > 2 ){
    G4cout << "PSDI Function Point A | 2D safety calculated to be: " << the2DSafety << G4endl;
  }
  
  //If we are within epsilon of a boundary, the next step should be made directly into the boundary so that boundary
  //processes can run. We will therefore not randomize the direction of the next step. Instead, we need to find the direction
  //toward the boundary in order to figure out how to angle the next vector. Note that an exception to this runs if we're *currently*
  //on a boundary (i.e. we're leaving from the boundary and the step is at a steep enough angle not to leave the epsilon region around
  //the surface.
  if( the2DSafety < epsilon ){

    //Debugging
    if( verboseLevel > 2 ){
      G4cout << "PSDI Function Point B | safety is smaller than epsilon, and finding direction to nearby boundary." << G4endl;      
    }    
    G4ThreeVector returnDir = FindDirectionToNearbyBoundary(track,the2DSafety);
    fParticleChange.ProposeMomentumDirection(returnDir.unit());
  }
  //3. If we're not within epsilon, then just purely randomize.
  else{

    //Debugging
    if( verboseLevel > 2 ){      
      G4cout << "PSDI Function Point B | safety is larger than epsilon, and we're just randomizing the next direction." << G4endl;      
    }
   
    //Randomize the final state momentum
    G4ThreeVector returnDir = G4RandomDirection();
    returnDir.setZ(0);
    fParticleChange.ProposeMomentumDirection(returnDir.unit());
  }
    
  ClearNumberOfInteractionLengthLeft();		// All processes should do this! 
  return &fParticleChange;

}

//Using the information about the track and the safety computed at the track point, we can find a direction to the
//nearby boundary.
G4ThreeVector G4CMPBogoliubovQPRandomWalkTransport::FindDirectionToNearbyBoundary(const G4Track& track, const G4double the2DSafety ){

  //Debugging
  if( verboseLevel > 2 ){
    G4cout << "---------- G4CMPBogoliubovQPRandomWalkTransport::FindDirectionToNearbyBoundary ----------" << G4endl;
  }
  
  //First, find the momentum direction. Since in the last step we went from outside epsilon to inside epsilon, this track
  //vector should be pointed roughly toward the boundary, rather than away from it. Use the momentum direction to find another
  //point that's a smidge farther from the boundary than the current track point. Here we assume that this point is still within
  //the original touchable. If that's not true, things will break here. But a delta of 1nm should be fine for this.
  G4ThreeVector momDir = track.GetMomentumDirection();

  //Sanity checks:
  //1. We set a default distance for computing another step and using this to find the direction to the wall as 1 nm
  G4double deltaPath = 1*nm;

  //2. However, if the2DSafety is less than 1nm, we will reset this to be the 2D safety so that we don't have weird edge
  //   cases where the upcoming shift sends a point across the surface.
  if( deltaPath > the2DSafety ) deltaPath = the2DSafety*0.999;

  //3. Now we shift the point by a bit and re-find the safety.
  G4ThreeVector shiftedPoint = track.GetPosition() - deltaPath*momDir;
  G4double shiftedPoint2DSafety = G4CMP::Get2DSafety(track.GetStep()->GetPreStepPoint()->GetTouchable(),
						     shiftedPoint,
						     track.GetMomentumDirection(),
						     false);
  //Debugging
  if( verboseLevel > 2 ){
    G4cout << "FDTNB Function Point A | pos: " << track.GetPosition() << ", shiftedPoint: " << shiftedPoint << ", original safety: " << the2DSafety << ", shiftedPoint2Dsafety: " << shiftedPoint2DSafety << G4endl;
  }

  //We now have two points: track.GetPosition() and shiftedPoint, and two safeties: the2DSafety and shiftedPoint2DSafety
  //We can now use these points to identify an angle between the current momDir and the surface normal. For now we'll assume
  //that this angle is in XY, but later (REL) we should come back and fix this to be more plane-agnostic.
  G4double deltaDistToSurface = shiftedPoint2DSafety - the2DSafety;

  //If the deltaDistToSurface is negative, this can arise in a few ways:
  //1. Complicated geometries with corners, though this is unlikely given how small the deltaPath parameter is -- would be very unlucky
  //   to get this.
  //2. If the step started on a boundary and launch our track at a steep angle, and we never *leave* the nearby epsilon region.
  //   Here, to re-launch ourselves back to the boundary, the final direction should be -1*momDir, with the attempted rotations
  //   down below. So if deltaDistToSurface is negative, we reverse the direction of momDir.
  
  //Sanity check
  if( deltaDistToSurface < 0 ){
    momDir = -1*momDir;
    G4ExceptionDescription msg;
    msg << "The DeltaDistToSurface is negative for some reason? This is probably an edge case that is actually ok, but flagging now for debugging.";
    G4Exception("G4CMPBogoliubovQPRandomWalkTransport::FindDirectionToNearbyBoundary", "BogoliubovQPRandomWalkTransport00X",JustWarning, msg);
  }


  G4double theta = acos(fabs(deltaDistToSurface)/deltaPath);

  //Debugging
  if( verboseLevel > 2 ){
    G4cout << "FDTNB Function Point B | theta for rotation: " << theta << G4endl;
  }

  //Now, we have the theta, but this theta could be in either direction relative to the surface. So what we'll do is try rotating our momDir vector
  //by theta in both directions, stepping the point forward by each, and re-finding safeties. This calculation is specific to the XY plane and
  //should be generalized (REL).
  G4ThreeVector option1 = momDir.rotateZ(theta);
  G4ThreeVector option2 = momDir.rotateZ(-2*theta);
  momDir.rotateZ(theta);

  //Debugging
  if( verboseLevel > 2 ){
    G4cout << "FDTNB Function Point C | potential direction-toward-boundary 1 (option1): " << option1 << G4endl;
    G4cout << "FDTNB Function Point C | potential direction-toward-boundary 2 (option2): " << option2 << G4endl;
  }

  //Recheck safety. Here, we compute the smaller of the two original safeties and make sure we set our
  //"probe" vectors' magnitudes no larger than that. (Here, the 0.9 is to keep it from being exactly one of the safeties)
  G4double smallerSafety = the2DSafety;
  if( shiftedPoint2DSafety < the2DSafety ) smallerSafety = shiftedPoint2DSafety;
  smallerSafety *= 0.9; 
  G4ThreeVector newPosOption1 = track.GetPosition() + smallerSafety*option1;
  G4ThreeVector newPosOption2 = track.GetPosition() + smallerSafety*option2;

  //Debugging
  if( verboseLevel > 2 ){
    G4cout << "FDTNB Function Point D | new probe position option 1: " << newPosOption1 << G4endl;
    G4cout << "FDTNB Function Point D | new probe position option 2: " << newPosOption2 << G4endl;
  }

  //Now compute the safeties of our two probe vectors' endpoints to see which is closer.
  G4double option1Safety = G4CMP::Get2DSafety(track.GetStep()->GetPreStepPoint()->GetTouchable(),
					      newPosOption1,
					      track.GetMomentumDirection(),
					      false);
  G4double option2Safety = G4CMP::Get2DSafety(track.GetStep()->GetPreStepPoint()->GetTouchable(),
					      newPosOption2,
					      track.GetMomentumDirection(),
					      false);

  //Debugging
  if( verboseLevel > 2 ){    
    G4cout << "FDTNB Function Point E | new probe position option 1 safety is " << option1Safety << G4endl;
    G4cout << "FDTNB Function Point E | new probe position option 2 safety is: " << option2Safety << G4endl;
  }

  //If option 1 safety is lower, it means we return option 1 as the direction to the boundary
  if( option1Safety < option2Safety ){
    if( verboseLevel > 2 ){
      G4cout << "FDNB Function Point F | Direction to the boundary is option 1: " << option1 << G4endl;
    }
    return option1;
  }
  else if( option2Safety < option1Safety ){
    if( verboseLevel > 2 ){
      G4cout << "FDNB Function Point G | Direction to the boundary is option 2: " << option2 << G4endl;
    }
    return option2;
  }
  else{
    G4ExceptionDescription msg;
    msg << "both option1 and option2 give the same distance for some reason? This seems to be an edge case.";
    G4Exception("G4CMPBogoliubovQPRandomWalkTransport::FindDirectionToNearbyBoundary", "BogoliubovQPRandomWalkTransport00X",FatalException, msg);
    G4ThreeVector dummy(0,0,0);
    return dummy;
  }
}
					     



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4CMPBogoliubovQPRandomWalkTransport::GetContinuousStepLimit(
                                       const G4Track& track,
                                       G4double previousStepSize,
                                       G4double currentMinimalStep,
                                       G4double& currentSafety)
{
  //Debugging
  if( verboseLevel > 2 ){
    G4cout << "---------- G4CMPBogoliubovQPRandomWalkTransport::GetContinuousStepLimit ----------" << G4endl;  
    G4cout << "GCSL Function Point A | In GetContinuousStepLimit." << G4endl;
  }
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
              const G4Track& track, G4double previousStepSize, G4ForceCondition* condition)
{
  //Debugging
  if( verboseLevel > 2 ){
    G4cout << "---------- G4CMPBogoliubovQPRandomWalkTransport::GetMeanFreePath ----------" << G4endl;
    G4cout << "GMFP Function Point A | track volume: " << track.GetVolume()->GetName() << G4endl;
    G4cout << "GMFP Function Point A | previousStepSize: " << previousStepSize << G4endl;
    G4cout << "GMFP Function Point A | momentum direction: " << track.GetMomentumDirection() << G4endl;
    G4cout << "GMFP Function Point A | position: " << track.GetPosition() << G4endl;
    G4cout << "GMFP Function Point A | global time: " << track.GetGlobalTime() << G4endl;
  }
  
  //This needs to be done so that we can update the SCUtils information, and since GetMeanFreePath for the discrete bit of this process should (?)
  //run first, it will hopefully happen before the AlongStepGPIL runs and requests the gap.
  if (UpdateMeanFreePathForLatticeChangeover(track)){
    UpdateSCAfterLatticeChange();
  }

  //Set up safety helper
  if(!fSafetyHelper) {
    G4TransportationManager* transportMgr ;
    transportMgr = G4TransportationManager::GetTransportationManager() ;
    fSafetyHelper = transportMgr->GetSafetyHelper();        
    fSafetyHelper->InitialiseHelper();
  }
  
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
  G4double energy = track.GetKineticEnergy();
  G4double velocity = track.GetVelocity();

  //Debugging
  if( verboseLevel > 2 ){
    G4cout << "GMFP Function Point B | CurrentVolume by track: " << currentVolume->GetName() <<G4endl;
    G4cout << "GMFP Function Point B | CurrentVolumePlusEps, by GetVolumeAtPoint: " << currentVolPlusEps->GetName() << G4endl;
    G4cout << "GMFP Function Point B | CurrentVolumeMinusEps, by GetVolumeAtPoint: " << currentVolMinEps->GetName() << G4endl;
    G4cout << "GMFP Function Point B | StepStatus: " << theStatus << G4endl;
  }

  //Initialize the lattice manager
  G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();
  
  //If the volume plus epsilon (in the direction of the step) and the current volume are not the same (and if we're
  //in a boundary-limited step), then we're in a turnaround step. This assumes that the real distance to the next
  //geometric feature is not smaller than the boundTolerance -- in that case, this may break down a bit. Return
  //that this is inactive.
  fTrackOnBoundary = false;
  if( currentVolPlusEps != currentVolume && theStatus == fGeomBoundary ){

    //NOTE: this above logic may run into issues in internal corners, if the next direction isn't pointed back into the volume. I.e. if a nm
    //pushes us across a corner back into World, we'll have an issue.
    
    //Debugging
    if( verboseLevel > 2 ){      
      G4cout << "GMFP Function Point C | In a turnaround step. Killing the transport GPIL." << G4endl;
    }    
    fTrackOnBoundary = true;
    isActive = false;
    return DBL_MAX;
  }

  //We can be on a boundary and still be "active" if we're in the volume into which we're moving.
  if( theStatus == fGeomBoundary ){
    
    //Debugging
    if( verboseLevel > 2 ){      
      G4cout << "GMFP Function Point D | On a boundary but not in a turnaround step." << G4endl;
    }
    fTrackOnBoundary = true;
  }

  //Debugging
  if( verboseLevel > 2 ){
    G4cout << "GMFP Function Point E | Not in a turnaround step. Continuing the transport GPIL." << G4endl;
    G4cout << "GMFP Function Point E | theStatus: " << theStatus << G4endl;
    G4cout << "GMFP Function Point E | Gap energy (drawn from SCUtils): " << fGapEnergy << G4endl;
    G4cout << "GMFP Function Point E | Energy (drawn from SCUtils): " << energy << G4endl;
    G4cout << "GMFP Function Point E | Dn (drawn from SCUtils): " << fDn << G4endl;
    G4cout << "GMFP Function Point E | Teff (drawn from SCUtils): " << fTeff << G4endl;
  }

  //Verify that, if we're on a boundary, the direction of momentum is going into a volume with a lattice. Since we've updated the
  //lattice at this point, we should confirm that the current procUtils lattice is not null, and that it is the
  //same as the lattice corresponding to what the latticeManager sees as belonging to currentVolPlusEps.
  //This is purely a check.
  if( theStatus == fGeomBoundary && !LM->HasLattice(currentVolPlusEps) ){ 
    G4ExceptionDescription msg;
    msg << "We're on a boundary (i.e. our boundary is now behind us) and find no lattice in the region where a BogoliubovQP is intending to go during the QP random walk transport's GetMeanFreePath function.";
    G4Exception("G4CMPBogoliubovQPRandomWalkTransport::GetMFP", "BogoliubovQPRandomWalkTransport001",FatalException, msg);
  }
  
  //Compute a step length corresponding to the nearest surface in 2D
  if ((energy>=fGapEnergy) && (fDn)>0){ 
    isActive = true;

    //Need to split this up into two scenarios: one where we're in bulk and one where we're on a boundary
    //Boundary case: here, we're not in a turnaround step but in the following boundary step where we're in
    //the volume that we'll continue into
    G4double the2DSafety;
    if( fTrackOnBoundary ){
      the2DSafety = G4CMP::Get2DSafety(track.GetStep()->GetPreStepPoint()->GetTouchable(),
				       track.GetPosition(),
				       track.GetMomentumDirection(),
				       true);
    }
    //Bulk case: simpler
    else{
      the2DSafety = G4CMP::Get2DSafety(track.GetStep()->GetPreStepPoint()->GetTouchable(),
				       track.GetPosition(),
				       track.GetMomentumDirection(),
				       false);
    }
    f2DSafety = the2DSafety;
    
    //Calculate the energy dependent diffusion constant
    G4double E_ratio = fGapEnergy/energy;
    G4double E_ratio2 = pow(E_ratio,2.0);
    fDiffConst = fDn*sqrt(1-E_ratio2);

    //Debugging
    if( verboseLevel > 2 ){      
      G4cout << "GMFP Function Point F | Diffusion constant (energy-adjusted): " << fDiffConst << G4endl;
    }
      
    
    //Using this diffusion constant, draw a time to reach the boundary
    //established by the safety. For now we use something simple and wrong to test,
    //taking the 2D safety to be the sigma of the diffusion broadening
    //but need the true expression from Grebenkov
    G4double timeStepToBoundary = SampleTimeStepFromFirstPassageDistribution(the2DSafety);    
    fTimeStepToBoundary = timeStepToBoundary;

    //Debugging
    if( verboseLevel > 2 ){      
      G4cout << "GMFP Function Point G | Time step to boundary = " << fTimeStepToBoundary << G4endl;
    }
    
    //With the time step to the boundary, we need to compute a "real," i.e. "full" distance
    //that can be compared to the other discrete processes. Since those use a characteristic time
    //multipled by a velocity that is quite high, we need to do the same.
    G4double thisMFP = timeStepToBoundary * velocity;    

    //There are some edge cases, where we have very small distances to the boundary due to vectors that "almost
    //pointed directly at the boundary." In these edge cases, (which may be resolved when we get the time of
    //first passage calculations in), the MFP as calculated with the above line will be less than the 2D safety
    //used as input. This is because of the quadratic form of the calculated timeStepToBoundary. In this scenario,
    //just manually set the MFP equal to the safety. (Usually, the MFP is much larger than the 2D safety.)
    //REL want to check that this is only for small enough distances not to matter, eventually.
    fVerySmallStep = false;
    if( thisMFP < the2DSafety ){
      fVerySmallStep = true;
      thisMFP = the2DSafety*fBoundaryFudgeFactor; //Fudge factor needs to go here as well because AlongStep GPIL uses PhysicalStep, which comes from PostStepGPIL
    }
    
    //Now, as a last measure, we set the number of interaction lengths left for this process to 1
    //and return a mean free path equal to this MFP, so that the only way that another process will
    //win over this in the discrete GPIL race is if the time to that process is less than the time
    //to this one.
    theNumberOfInteractionLengthLeft = 1;

    //Debugging
    if( verboseLevel > 2 ){      
      G4cout << "GMFP Function Point H | Setting the number of interaction lengths for RWTransport's Discrete component to 1, and this MFP = " << thisMFP << G4endl;
    }
    return thisMFP;
  }
  else{
    //This indicates something aphysical.
    G4ExceptionDescription msg;
    msg << "QP energy is too low or we're missing a Dn. Returning DBL_MAX for GetMFP.";
    G4Exception("G4CMPBogoliubovQPRandomWalkTransport::GetMeanFreePath", "BogoliubovQPRandomWalkTransport002",FatalException, msg);
    isActive = false;
    return DBL_MAX;
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//This is an amended version of the logic block in G4CMPVProcess that allows us to update the lattice in procUtils and in SCUtils when it changes.
//Since we don't use a rateModel or other things explicitly here, it's not the exact same code, which is why I'm okayish with havine a separate
//funciton in this class.
G4bool G4CMPBogoliubovQPRandomWalkTransport::UpdateMeanFreePathForLatticeChangeover(const G4Track& aTrack)
{
  //Debugging
  if( verboseLevel > 2 ){
    G4cout << "---------- G4CMPBogoliubovQPRandomWalkTransport::UpdateMeanFreePathForLatticeChangeover ----------" << G4endl;
    G4cout << "UMFPFLC Function Point A | Loading data for track after lattice changeover, process: " << this->GetProcessName() << G4endl;
    G4cout << "UMFPFLC Function Point A | Track length: " << aTrack.GetTrackLength() << G4endl;
    G4cout << "UMFPFLC Function Point A | Current lattice a la lattice manager: " << G4LatticeManager::GetLatticeManager()->GetLattice(aTrack.GetVolume()) << G4endl;
  }
    
  //Always do a check to see if the current lattice stored in this process is equal to the one that represents
  //the volume that we're in. Note that we can't do this with the "GetLattice()" and "GetNextLattice()" calls
  //here because at this point in the step, the pre- and post-step points both point to the same volume. Since
  //GetMeanFreePath is run at the beginning, I think the point at which a boundary interaction is assessed comes
  //later (hence why we can use that info in PostStepDoIts but not here.) Adding a statement about track length here,
  //since it seems that when a particle spawns it doesn't necessarily trigger this block, and I think we want it to.
  if( (((this->theLattice) && G4LatticeManager::GetLatticeManager()->GetLattice(aTrack.GetVolume())) &&
       (this->theLattice != G4LatticeManager::GetLatticeManager()->GetLattice(aTrack.GetVolume()))) ||
      aTrack.GetTrackLength() == 0.0 ){

    //Debugging
    if( verboseLevel > 2 ){
      G4cout << "UMFPFLC Function Point B | Step length associated with this is " << aTrack.GetStep()->GetStepLength() << G4endl;
      G4cout << "UMFPFLC Function Point B | Successfully changed over to a new lattice for process " << this->GetProcessName() << G4endl;
    }
        
    //REL noting that if physical lattices are not 1:1 with volumes, something may get broken here... Should check a scenario of segmented SC...    
    this->LoadDataForTrack(&aTrack);
    return true;
    
  }

  //Debugging
  if( verboseLevel > 2 ){
    G4cout << "UMFPFLC Function Point C | Did not successfully change over to a new lattice for process " << this->GetProcessName() << G4endl; 
  }
  return false;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//This is meant to update superconductor info for the process if we move into a new lattice. Custom to this process, which does not have an associated
//rate model.
void G4CMPBogoliubovQPRandomWalkTransport::UpdateSCAfterLatticeChange()
{
  //Debugging
  if( verboseLevel > 2 ){
    G4cout << "---------- G4CMPBogoliubovQPRandomWalkTransport::UpdateSCAfterLatticeChange ----------" << G4endl;
    G4cout << "USCALC Function Point A | Updating SC After Lattice Change" << G4endl;
  }
  
  //First, determine if the new lattice is a SC. If not, then set the SCUtils info to null for this process  
  if( (this->theLattice)->GetSCDelta0() <= 0 ){
    this->SetCurrentSCInfoToNull();
    return;
  }

  //If it is a SC, then we should update the SC information for the SC utils class within the base of this.
  //Also, handle the checking/updating of the lookup tables to be used for each SC.
  this->LoadLatticeInfoIntoSCUtils(this->theLattice);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//Sample a time step for first passage through a 2D disk of radius the2DSafety. This is from the 2013 Grebenkov paper on
//efficient monte carlo techniques for diffusion. (arXiv:1304.7807, physics.comp-ph)
G4double G4CMPBogoliubovQPRandomWalkTransport::SampleTimeStepFromFirstPassageDistribution(G4double the2DSafety)
{
  //Hardcoded -- can change if we want to
  std::string samplingTechnique = "AcceptanceRejection";
  
  //Use an acceptance-rejection technique with the (known) pdf
  if( samplingTechnique == "AcceptanceRejection" ){
    G4double dimensionlessTime = SampleDimensionlessTimeStepUsingAcceptanceRejectionTechnique();
    return dimensionlessTime * the2DSafety*the2DSafety / fDiffConst;
  }

  //Spot for other techniques that are more efficient (Inversion/Newton, etc.)
  else{
    G4ExceptionDescription msg;
    msg << "Attempting to use an unknown time sampling technique, " << samplingTechnique << " for first-passage time. Throwing an error.";
    G4Exception("G4CMPBogoliubovQPRandomWalkTransport::SampleTimeeStepFromFirstPassageDistribution", "BogoliubovQPRandomWalkTransport00X",FatalException, msg);
  }
}

//This samples a dimensionless time step using A/R technique. 
G4double G4CMPBogoliubovQPRandomWalkTransport::SampleDimensionlessTimeStepUsingAcceptanceRejectionTechnique()
{
  if( verboseLevel > 2 ){
    G4cout << "---------- G4CMPBogoliubovQPRandomWalkTransport::SampleDimensionlessTimeStepUsingAcceptanceRejectionTechnique ----------" << G4endl;
    G4cout << "SDTSUART Function Point A | Sampling dimensionless time step with A/R technique." << G4endl;
  }
  
  //Define min and max allowable sampled dimensionless times. These are
  //set somewhat arbitrarily but are basically meant to ensure that we
  //aren't including in our domain regions where the pdf is vanishingly small
  //(which increases our execution time).
  G4double minT = 0.001;
  G4double maxT = 3;
  G4double maxVal = 4.06615; //Hardcoded REL, but the max value of our calculated Grebenkov disk pdf for dSc/dt

  //Loop indefinitely -- when we find a point under the curve, we'll break
  G4double sampledT = -1;
  while(1){
    //Sample a point in the 2D space spanned by this pdf
    double sampleT = G4UniformRand()*(maxT-minT);    
    double sampleY = G4UniformRand()*maxVal;

    //Compute the pdf at the given sample T and compare to the sampleY
    //These need to be made into less-hardcoded numbers.
    double BesselJ1_evalAt_alpha00 = 0.519153; //Hardcoded REL
    double alpha00 = 2.4048; //Hardcoded REL    
    double derivHigh = -2.0*alpha00/BesselJ1_evalAt_alpha00*exp(-1*alpha00*alpha00*sampleT);
    double derivLow = exp(-1.0/4.0/sampleT) * (-0.5/sampleT/sampleT + 0.5/sampleT - 16*sampleT);

    //Debugging
    if( verboseLevel > 2 ){
      G4cout << "SDTSUART Function Point B | Sampled T: " << sampleT << ", sampled Y: " << sampleY << ", -derivLow: " << -1*derivLow << ", -derivHigh: " << -1*derivHigh << G4endl;
    }

    //This defines the transition between the low-time approximation and the high-time approximation
    if( sampleT <= 0.15 ){
      if( sampleY < -1*derivLow ){
	sampledT = sampleT;
	break;
      }
    }
    else{
      if( sampleY < -1*derivHigh ){
	sampledT = sampleT;
	break;
      }
    }   
  }

  //Sanity check
  if( sampledT < 0 ){
    G4ExceptionDescription msg;
    msg << "Somehow failed to sample an appropriate time using the acceptance/rejection technique. Throwing an error.";
    G4Exception("G4CMPBogoliubovQPRandomWalkTransport::SampleDimensionlessTimeStepUsingAcceptanceRejectionTechnique", "BogoliubovQPRandomWalkTransport00X",FatalException, msg);
  }

  //fOutfile << sampledT << G4endl;  

  return sampledT;
  
}

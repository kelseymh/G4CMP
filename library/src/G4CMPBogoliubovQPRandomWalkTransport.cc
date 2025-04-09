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

  SetProcessSubType(fBogoliubovQPRandomWalkTransport);
  
  //This time step will be overwritten by the step-limiting length (discrete process GPIL race)
  fTimeStep = 0;
  fPathLength =  0.0;
  fPreDiffusionPathLength = 0.0;
  fDiffConst =  0.0;
  fBoundaryFudgeFactor = 1.0001;
  //fBoundaryFudgeFactor = 1.002;
  fEpsilonForWalkOnSpheres = 1*CLHEP::um;
  
  //fSafetyHelper is initialized in AlongStepGPIL
  fSafetyHelper=nullptr;
  
  //Temporary REL
  //fOutfile.open("/Users/ryanlinehan/QSC/Sims/Geant4/scRebuild-build/RandomWalkSampledDimlessTimes.txt",std::ios::trunc);

  fBoundaryHistoryTrackID = -1;
  fBoundaryHistory.clear();
  fMaxBoundaryHistoryEntries = 6; //HARDCODED REL -- FIX
  fQPIsStuck = false;
  fOutgoingSurfaceTangent1 = G4ThreeVector(0,0,0);
  fOutgoingSurfaceTangent2 = G4ThreeVector(0,0,0);
  fDotProductDefiningUniqueNorms = 0.99;
  fStuckInCornerThreshold = fEpsilonForWalkOnSpheres;
  fNeedSweptSafetyInGetMFP = false;
  
  //REL A GENTLE REMINDER THAT WE CANNOT INSTANTIATE CONFIG MANAGER VALUES ONLY IN PROCESS OR RATE CONSTRUCTORS BECAUSE THEY COME BEFORE
  //THE CONFIG MANAGER IS INSTANTIATED
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
  //Need to do this here because if we do it in the constructor, that runs before the ConfigManager instance
  //is defined in the current structure of the examples
  verboseLevel = G4CMPConfigManager::GetVerboseLevel();  
  
  //Debugging
  if( verboseLevel > 5 ){
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
  if(verboseLevel>5){
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
  if( verboseLevel>5){
    G4cout << "ASGPIL Function Point B | fPathLength = fPreDiffusionPathLength = currentMinimalStep: " << fPathLength << G4endl;
    G4cout << "ASGPIL Function Point B | velocity: " << track.GetVelocity() << G4endl;
  }

  //If we're in a turnaround step, then kill
  if( isActive == false ){

    //Some debugging
    if( verboseLevel>5 ){
      G4cout << "ASGPIL Function Point C | In a turnaround step. Killing the transport GPIL." << G4endl;
    }
    return DBL_MAX;
  }
  
  //Get some track and SC properties. The diffusion constant and other superconductor-specific information are accessed via the
  //SCUtils class that this inherits from. No need to actually grab it from somewhere: it's already a data member
  G4double energy = track.GetKineticEnergy();
  G4double velocity = track.GetVelocity();
  G4ThreeVector momentumDir = track.GetMomentumDirection();  
  if( verboseLevel > 5 ){
    G4cout << "ASGPIL Function Point D | Gap energy (drawn from SCUtils): " << fGapEnergy << ", particle energy: " << energy << G4endl;
    G4cout << "ASGPIL Function Point D | Dn (drawn from SCUtils): " << fDn << ", and fDiffConst: " << fDiffConst << G4endl;
    G4cout << "ASGPIL Function Point D | Teff (drawn from SCUtils): " << fTeff << G4endl;
  }
    
  //If our energy is appropriate and we can see diffusion info, trigger the "meat" of this function
  if ((energy>=fGapEnergy) && (fDn)>0){ 
    isActive = true;
    *selection = NotCandidateForSelection;

    //Debugging
    if( verboseLevel > 5 ){      
      G4cout << "ASGPIL Function Point E | isActive is true. currentMinimalStep/velocity: " << currentMinimalStep/velocity << ", fTimeStepToBoundary: " << fTimeStepToBoundary << G4endl;
    }

    //Here, we need to return a diffusion-displacement ("diffusion-folded") path length. If the currentMinimalStep matches
    //the fTimeStepToBoundary, it means that the winning discrete GPIL process is THIS process's discrete bit. In this case,
    //the diffusion displacement distance is simple: it's just the distance to the boundary (up to the fudge factor that we use
    //to actually trigger a boundary interaction in G4), and the time is just the time step to the boundary.
    double timeTolerance = 1E-10; //For floating point errors    
    if( (fabs(currentMinimalStep/velocity - fTimeStepToBoundary) < timeTolerance) || fVerySmallStep ){

      //Let's first handle the "very small step" scenarios.
      if( fVerySmallStep ){
	if( verboseLevel > 5 ){
	  G4cout << "ASGPIL Function Point F_0 | Looks like we're in a very small step. We are going to launch the track a bit farther than the 2D safety, and hopefully this gets us to a boundary. (Conditions for this include nearly-stuck QPs)" << G4endl;
	}
	fPathLength = f2DSafety * fBoundaryFudgeFactor;
	fTimeStep = fTimeStepToBoundary;
      }
      else{
       
	//If we're not very close to the boundary, then we just make our diffusion-folded path length equal to juuuuuust under the
	//distance to the boundary. That way transportation will never win this alongStepGPIL.
	//This also should trigger if we're starting from the boundary
	if( f2DSafety >= fEpsilonForWalkOnSpheres || (fTrackOnBoundary) ){
	  fPathLength = f2DSafety / fBoundaryFudgeFactor; 
	  fTimeStep = fTimeStepToBoundary;
	  
	  //Debugging
	  if( verboseLevel > 5 ){
	    G4cout << "ASGPIL Function Point F_1 | Looks like the boundary-limited case applies here, with 2DSafety >= epsilon. Returning fPathLength just under f2DSafety = " << fPathLength << G4endl;
	  }
	}
	
	//Otherwise, if we are close to the boundary, then in this step we need to make our diffusion-folded path length equal to juuuuuuust
	//larger than the distance to the boundary. This ensures we will always have transportation win and take us to the boundary.
	else{
	  fPathLength = f2DSafety * fBoundaryFudgeFactor;
	  fTimeStep = fTimeStepToBoundary;
	  
	  //Debugging
	  if( verboseLevel > 5 ){
	    G4cout << "ASGPIL Function Point F_2 | Looks like the boundary-limited case applies here, with 2DSafety < epsilon. Returning fPathLength just over f2DSafety = " << fPathLength << G4endl;
	  }
	}
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
      while( fPathLength >= f2DSafety );

      //Debugging
      if( verboseLevel > 5 ){
	G4cout << "ASGPIL Function Point G | Sigma1D: " << sigma1D << " from fDiffConst: " << fDiffConst << G4endl;
	G4cout << "ASGPIL Function Point G | Looks like a different discrete process wins the GPIL race. Returning fPathLength = " << fPathLength << G4endl;
      }
    }

    //Debugging
    if( verboseLevel > 5 ){
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
  if( verboseLevel > 5 ){
    G4cout << "---------- G4CMPBogoliubovQPRandomWalkTransport::PostStepGetPhysicalInteractionLength ----------" << G4endl;
    G4cout << "PSGPIL Function Point A | In PostStepGetPhysicalInteractionLength" << G4endl;
  }

  //Since we're overriding this function as well, we'll have to call the MFP one (since it's called here in the base class).
  //This is all just to get the SCUtils updated at the beginning of the step calculus. Doing it in MFP because that's where
  //it's "usually" done with all of the other classes (even though it doesn't really matter which step it's called in here)
  double mfp = GetMeanFreePath(track, previousStepSize, condition);
  fNeedSweptSafetyInGetMFP = false; //Reset this here so that if PostStepDoIt doesn't run, we can still reset.
  return mfp;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// "Stretch" version
G4VParticleChange* G4CMPBogoliubovQPRandomWalkTransport::AlongStepDoIt(const G4Track& track, const G4Step& step) {
  
  //Debugging
  if( verboseLevel > 5 ){
    G4cout << "---------- G4CMPBogoliubovQPRandomWalkTransport::AlongStepDoIt ----------" << G4endl;
    G4cout << "ASDI Function Point A | Step Status from track pre-step point: " << track.GetStep()->GetPreStepPoint()->GetStepStatus() << G4endl;
    G4cout << "ASDI Function Point A | Step Status from step pre-step point: " << step.GetPreStepPoint()->GetStepStatus() << G4endl;
    G4cout << "ASDI Function Point A | Step Status from step post-step point: " << step.GetPostStepPoint()->GetStepStatus() << G4endl;
  }

  //Particle change initialization and setting
  fParticleChange.Initialize(track);  
  fParticleChange.ProposeMomentumDirection(step.GetPostStepPoint()->GetMomentumDirection());
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
  if( verboseLevel > 5 ){
    G4cout << "ASDI Function Point B | stepStartGlobalTime: " << stepStartGlobalTime << G4endl;
    G4cout << "ASDI Function Point B | stepEndGlobalTime  : " << stepEndGlobalTime << G4endl;
    G4cout << "ASDI Function Point B | stepTransportOnlyDeltaT: " << stepTransportOnlyDeltaT << G4endl;
  }
  
  //As a reminder, the velocity for QPs is always a single number that stays constant.
  G4double velocity = step.GetPostStepPoint()->GetVelocity();
  fParticleChange.ProposeVelocity(velocity);
  fPositionChanged = false;
  G4double stepLength = step.GetStepLength();
  
  // Check if the particle met conditions to do random walk from GPIL command. Here, this occurs
  // during a turnaround step where we want to set step lengths to zero. If we don't set the proposedTrueStepLength
  // to zero, then we'll end up mucking up our boundary processes which expect a zero step length
  // during turnaround steps.
  if(!isActive) {
    fPathLength = stepLength;
    fPreDiffusionPathLength = stepLength;

    //Debugging
    if( verboseLevel > 5 ){
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
    if( verboseLevel > 5 ){
      G4cout << "ASDI Function Point D | velocity is set and proposed to " << velocity << ", which should never change." << G4endl;
    }
    fParticleChange.ProposeVelocity(velocity);

    
    //Have to stratify by three cases: off-boundary, on-boundary, and on-boundary-and-stuck
    //Case 1: off-boundary
    // - Here we just generate a purely random vector in 2D and move it by the fPathLength
    fOldPosition = step.GetPreStepPoint()->GetPosition();
    if( !fTrackOnBoundary ){

      //Sub case 1: farther than epsilon away from boundary? Yay -- get to randomize
      if( f2DSafety >= fEpsilonForWalkOnSpheres ){
	G4ThreeVector thisRandomDir = G4RandomDirection();
	thisRandomDir.setZ(0);
	fNewDirection = thisRandomDir.unit();

	//Debugging
	if( verboseLevel > 5 ){	
	  G4cout << "ASDI Function Point D_1 | Since we're in the bulk and >eps from a boundary, we're randomizing new direction to be " << fNewDirection << G4endl;
	}

	
      }
      //Otherwise, don't interfere with the direction needed to get the step to the boundary. Just
      //let the new direction be the old direction so we can respect the short hop to the boundary.
      else{
	fNewDirection = step.GetPreStepPoint()->GetMomentumDirection();

	//Debugging
	if( verboseLevel > 5 ){	
	  G4cout << "ASDI Function Point D_2 | Since we're in the bulk and <eps from a boundary, we're keeping the new direction as " << fNewDirection << G4endl;
	}
      }
    }
    else{    

      //Define a tolerance for the dot product of the emerging vector and the surface norm, and then get the surface norm (here, the momentum dir
      //if the boundary interactions are handled correctly.)

      //Want to generate our points such that we start from the norm and work outwards
      G4int nV = G4CMPConfigManager::GetSafetyNSweep2D();
      G4double dotProductThreshold_Norm = G4CMP::ComputeDotProductThreshold_Norm(nV);
      G4double dotProductThreshold_Tang = G4CMP::ComputeDotProductThreshold_Tang(nV);
      
      //G4double epsilonDotProductForNorm = 0.0896393089; //This matches the granularity of the from-boundary 2D safety. //0.07; //This ultimately needs to match the granularity of the 2DSafety depending on the curvature REL. In this
      //case we're talking about things that are being checked for being 86 degrees apart
      //G4double epsilonDotProductForTangent = 0.002; //Same, but in this case we're checking for things that are three degrees apart
      G4ThreeVector surfaceNorm = track.GetMomentumDirection();

      //Case 1: not stuck. In this case, we pull the pre-step point, whose vector should actually just be the surface normal in the
      //direction of the particle's new momentum (this is from the QP boundary class). We can then randomly generate a vector in 2D,
      //and accept that the distance is just that drawn from GetMFP
      if( !fQPIsStuck ){
	G4ThreeVector theNewDirection;
	do{
	  theNewDirection = G4RandomDirection();
	  theNewDirection.setZ(0);
	  theNewDirection = theNewDirection.unit();
	  //Sometimes there may be issues here if our particle is errantly on a top/bottom boundary
	}
	while( theNewDirection.dot(surfaceNorm) <= dotProductThreshold_Norm );
	fNewDirection = theNewDirection;
	
	//Debugging
	if( verboseLevel > 5 ){	
	  G4cout << "ASDI Function Point D_3 | Since we're on a surface, we're launching in the new direction " << fNewDirection << ", which is consistent with a dotProductThreshold_Norm of " << dotProductThreshold_Norm << G4endl;
	}
      }

      //Case 2: stuck. Here we have a different safety (and different distance) computed in GetMFP. We still need to modify the direction of the
      //track here, however.
      else{
	
	//With this safety, generate a new direction somewhere between the surface tangents. Probably can speed this up but for now
	//get something that works REL)
	G4ThreeVector theNewDir(0,0,0);
	G4double minDot = fOutgoingSurfaceTangent1.dot(fOutgoingSurfaceTangent2);
	G4double acosTang1Dot = 0;
	G4double acosTang2Dot = 0;
	G4double acosMinDot = 0;
	
	//Debugging
	if( verboseLevel > 5 ){	
	  G4cout << "ASDI Function Point D_4 | minDot is " << minDot << " between " << fOutgoingSurfaceTangent1 << " and " << fOutgoingSurfaceTangent2 << G4endl;
	  G4cout << "ASDI Function Point D_4 | theNewDir.dot(fOutgoingSurfaceTangent1): " << theNewDir.dot(fOutgoingSurfaceTangent1) << ", theNewDir.dot(fOutgoingSurfaceTangent2): " << theNewDir.dot(fOutgoingSurfaceTangent2) << G4endl;
	}

	
	do{
	  theNewDir = G4RandomDirection();
	  theNewDir.setZ(0);
	  theNewDir = theNewDir.unit();

	  //Calculate some angles
	  acosTang1Dot = acos(theNewDir.dot(fOutgoingSurfaceTangent1));
	  acosTang2Dot = acos(theNewDir.dot(fOutgoingSurfaceTangent2));
	  acosMinDot = acos(minDot);
	}

	//The logic here is a bit mucky, but basically all of the following conditions must be met:
	//1. The new dir should have a larger dot product with both of the tangent vectors than the two tangent vectors have with each other.
	//2. Since in the case of angles >120 degrees, area opens up on the "wrong" side of the corner that satisfies condition 1, we also need
	//   to satisfy the criterion that the sum of the angles between the new dir and each dir must be less than 180 degrees.
	//   This 180 degree absolute angular scale assumes that only at concave corners will we potentially get stuck, and for convex corners
	//   we won't ever arrive here. I can't imagine this is violated?
	//3. Lastly, also need to make sure that we're not overly close to either of the tangent vectors (don't want to launch off parallel to a
	//   surface for computational efficiency's sake).
	while( !(theNewDir.dot(fOutgoingSurfaceTangent1) > minDot && theNewDir.dot(fOutgoingSurfaceTangent2) > minDot &&
		 (acosTang1Dot + acosTang2Dot < CLHEP::pi) &&
		 theNewDir.dot(fOutgoingSurfaceTangent1) < dotProductThreshold_Tang && theNewDir.dot(fOutgoingSurfaceTangent2) < dotProductThreshold_Tang) );
	
	
	
	//Now set the direction
	fNewDirection = theNewDir;
      }
    }
    fNewPosition = fOldPosition + fPathLength*fNewDirection;

    /*
    //There are apparently some scenarios in which the returned position of the QP is actually on the top or bottom of the film? 
    //Anyway, just check that this isn't happening, and if it is, nudge the QP away from the film surface.
    G4double safetyInZ = G4CMP::GetSafetyInZ(track.GetStep()->GetPreStepPoint()->GetTouchable(),track.GetStep->GetPreStepPoint);
    if( safetyInZ )
    */

    
    
    //Debugging
    if( verboseLevel > 5 ){
      G4cout << "ASDI Function Point E | fOutgoingSurfaceTangent1: " << fOutgoingSurfaceTangent1 << G4endl;
      G4cout << "ASDI Function Point E | fOutgoingSurfaceTangent2: " << fOutgoingSurfaceTangent2 << G4endl;
      G4cout << "ASDI Function Point E | time of pre-step point is: " << step.GetPreStepPoint()->GetGlobalTime() << G4endl;
      G4cout << "ASDI Function Point E | time of post-step point is: " << step.GetPostStepPoint()->GetGlobalTime() << G4endl;
      G4cout << "ASDI Function Point E | fPathLength selected: " << fPathLength << G4endl;
      G4cout << "ASDI Function Point E | fNewDirection: " << fNewDirection << G4endl;
    }

    //Test. I think CheckNextStep may also be valuable here, especially when we get into scenarios where we enter daughter volumes
    double nextStepSafety = 0;    
    double nextStepLength = fSafetyHelper->CheckNextStep(fOldPosition,fNewDirection,fPathLength,nextStepSafety);

    //Debugging
    if( verboseLevel > 5 ){
      G4cout << "ASDI Function Point F | Checking next step. CheckNextStep's next step length: " << nextStepLength << ", nextStepSafety: " << nextStepSafety << G4endl;
    }

    //We're actually going to try using checkNextStep instead. Here, if we're not on a boundary, the returned step should be kInfinity and
    //the safety should be nonzero.
    if( nextStepLength == kInfinity ){ //We're out in the bulk

      //Debugging
      if( verboseLevel > 5 ){
	G4cout << "ASDI Function Point G | the proposed step does not cross a boundary. Setting position manually." << G4endl;
      }
      fSafetyHelper->ReLocateWithinVolume(fNewPosition);
      fParticleChange.ProposeMomentumDirection(fNewDirection); 

      //Since we are forced to use the pre-step point's velocity for propagation of this step (and the pre-step point's velocity is
      //the one we started with), transportation will add a time corresponding to traveling the calculated path length at that velocity.
      //We should subtract that time off our final proposed time. First, debugging.
      if( verboseLevel > 5 ){
	G4cout << "ASDI Function Point H | fTimeStep: " << fTimeStep << ", timeChangefromTransportationOnly: " << stepTransportOnlyDeltaT << G4endl;
	G4cout << "ASDI Function Point H | timeChangeFromTransportationOnly: " << stepTransportOnlyDeltaT << G4endl;
      }
      fParticleChange.ProposeLocalTime(fTimeStep-stepTransportOnlyDeltaT);


      //I think this needs to be set to the step-limiting *old, pre-diffusion* path length, since other processes that aren't the step-limiting one will
      //need to subtract off a distance. That distance basically needs to be velocity * deltaT. The current process that limits the step
      //(here, not transportation) will zero out its number of interaction lengths. First, debugging
      if( verboseLevel > 5 ){
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
	msg << "Somehow the CheckNextStep returned a step length that is not kInfinity but the step status thinks it's not fGeomBoundary. It seems we may have misjudged the distance to our boundary.";
	G4Exception("G4CMPBogoliubovQPRandomWalkTransport::AlongStepDoIt", "BogoliubovQPRandomWalkTransport00X",FatalException, msg);
      }
      else{

	//Debugging
	if( verboseLevel > 5 ){
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
  if( verboseLevel > 5 ){
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
  if( verboseLevel > 5 ){
    G4cout << "---------- G4CMPBogoliubovQPRandomWalkTransport::PostStepDoIt ----------" << G4endl;
  }

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
  
  //Determine if we end up close to a boundary using Get2DSafety. Here it takes the post-step point (track.GetPosition())
  G4double the2DSafety = G4CMP::Get2DSafety(track.GetStep()->GetPreStepPoint()->GetTouchable(),
					    track.GetPosition(),
					    track.GetMomentumDirection(),
					    false);


  //There are occasional edge cases where the 2D safety is zero here, which it really shouldn't be. I believe these are probably
  //because the safety checks that transportation uses to determine if it wins the alongStepDoIt race may try to be efficient
  //and use the simple DistToIn functions (and are not limited by kCarTolerance), while my directional-based checks to mother volumes
  //are limited by kCarTolerance and can be zero if they're below that threshold (even if they're not ACTUALLY zero).
  //(I think, that is -- hard to suss these errors out.) Here, we try shifting the particle backwards along its trajectory by
  //a small distance and just return that.
  if( the2DSafety == 0 ){
    G4cout << "In PSDI, Somehow our initial the2DSafty in PostStepDoIt is zero... " << G4endl;
    G4cout << "Track position: " << track.GetPosition() << ", momentum direction: " << track.GetMomentumDirection() << G4endl;
    
    
    //G4ExceptionDescription msg;
    //msg << "Killing.";
    //G4Exception("G4CMPBogoliubovQPRandomWalkTransport::PostStepDoIt", "BogoliubovQPRandomWalkTransport00X",FatalException, msg);

    //Get displacement (diffusion-folded) length of step and move ourselves back by an amount that depends on that length
    G4ThreeVector newPos;
    G4double stepLength = (track.GetStep()->GetPreStepPoint()->GetPosition()-track.GetStep()->GetPostStepPoint()->GetPosition()).mag();
    if( stepLength > fEpsilonForWalkOnSpheres ){
      newPos = track.GetPosition()-track.GetMomentumDirection()*fEpsilonForWalkOnSpheres;
    }
    else{
      newPos = track.GetPosition()-track.GetMomentumDirection()*0.5*stepLength;
    }
    G4ThreeVector rejiggeredMomDir = G4RandomDirection();
    rejiggeredMomDir.setZ(0);
    rejiggeredMomDir = rejiggeredMomDir.unit();
    
    double nextStepSafety = 0;
    double nextStepLength = fSafetyHelper->CheckNextStep(newPos,rejiggeredMomDir,DBL_MAX,nextStepSafety);
    fSafetyHelper->ReLocateWithinVolume(newPos);

    G4StepStatus theStatus = track.GetStep()->GetPostStepPoint()->GetStepStatus();
    
    G4cout << "In our 2dSafety=0 scenario, nextStepLength for newPos after adjustment is " << nextStepLength << " and the post-step status is: " << theStatus << G4endl;

    //G4ExceptionDescription msg;
    //msg << "Killing.";
    //G4Exception("G4CMPBogoliubovQPRandomWalkTransport::PostStepDoIt", "BogoliubovQPRandomWalkTransport00X",FatalException, msg);

    fParticleChange.ProposePosition(newPos);
    fParticleChange.ProposeMomentumDirection(rejiggeredMomDir);
    ClearNumberOfInteractionLengthLeft();		// All processes should do this! 
    return &fParticleChange;

  }  

  
  //Debugging
  if( verboseLevel > 5 ){
    G4cout << "PSDI Function Point A | 2D safety calculated to be: " << the2DSafety << G4endl;
  }
  
  //If we are within epsilon of a boundary, the next step should be made directly into the boundary so that boundary
  //processes can run. We will therefore not randomize the direction of the next step. Instead, we need to find the direction
  //toward the boundary in order to figure out how to angle the next vector. Note that an exception to this runs if we're *currently*
  //on a boundary (i.e. we're leaving from the boundary and the step is at a steep enough angle not to leave the epsilon region around
  //the surface. //REL This may not run if another process sends us within epsilon of the boundary... do we need to do anything about this?
  if( the2DSafety < fEpsilonForWalkOnSpheres ){

    //Debugging
    if( verboseLevel > 5 ){
      G4cout << "PSDI Function Point B | safety is smaller than epsilon, and finding direction to nearby boundary." << G4endl;      
    }
    G4bool needToRepeat = false;
    G4ThreeVector returnDir = FindDirectionToNearbyBoundary(track,the2DSafety,needToRepeat,false);

    //    G4cout << "Direction to nearby boundary: " << returnDir << G4endl;
    
    //Sometimes the basic G4 DistToIn functions for calculating safety to daughters don't return a mathematically correct distance to the
    //inside of the daughter volume (see G4Box, G4Trd DistToIn for examples). In this scenario, the FindDirectionToNearbyBoundary will throw
    //a flag and identify that something about the boundary calculation is funky. In this case, we will need to recompute the 2D safety and demand that if
    //the safety is to a daughter, we compute it using the sweep method (rather than the "fast" and sometimes wrong directionless  DistToIn).
    //We'll do this for the initial Get2DSafety (i.e. Get2DSafetyWithDirection), which should also return a direction. This is a bit of a nuclear option
    //since if there are lots of daughters, this will compute a sweep for AAAAAALLL of them. So we only want to run this sparingly.
    //Also happens if we're near corners where the findDirection strategy math breaks down
    if( needToRepeat ){
      G4cout << "Need to repeat tripped." << G4endl;
      std::pair<G4double,G4ThreeVector> the2DSafetyAndDir = G4CMP::Get2DSafetyWithDirection(track.GetStep()->GetPreStepPoint()->GetTouchable(),
											    track.GetPosition(),
											    track.GetMomentumDirection(),
											    false);
      the2DSafety = the2DSafetyAndDir.first;
      G4ThreeVector safetyDir = the2DSafetyAndDir.second;

      //Do a second check, moving along the direction of the safety
      G4ThreeVector shiftedPosForTest = track.GetPosition()+safetyDir*the2DSafety*0.5;
      std::pair<G4double,G4ThreeVector> the2DSafetyAndDir_Shifted = G4CMP::Get2DSafetyWithDirection(track.GetStep()->GetPreStepPoint()->GetTouchable(),
												    shiftedPosForTest,
												    track.GetMomentumDirection(),
												    false);

      //Another sanity check -- if we've repeated and STILL fail to pick a direction that gives an expected change in the safety when we
      //shift along the nominal safety direction, then something is wrong.
      if( fabs((the2DSafetyAndDir_Shifted.first / the2DSafety) - 0.5) > 0.1 ){

	//Edge case: if the original safety is not zero and the checked safety is zero, this direction is probably okay, but this check metric
	//will return -0.5. Circumvent this.
	if( the2DSafetyAndDir_Shifted.first == 0 && the2DSafety > 0 ){
	  //Flag
	  
	}
	else{
	  G4ExceptionDescription msg;
	  msg << "After first needToRepeat flagged, we are checking the repeated safety/dir calculation and find that the safetyDir is not plausibly in the direction of the nearest boundary.";
	  G4Exception("G4CMPBogoliubovQPRandomWalkTransport::PostStepDoIt", "BogoliubovQPRandomWalkTransport00X",FatalException, msg);
	}
      }


      
      //In some scenarios, the original 2D safety is below fEpsilonForWalkOnSpheres and the new, correct 2D safety is not. If this is the case, we will
      //likely not actually hit the surface after pointing ourselves toward it (since we make our step length only slightly larger than the
      //original computed safety. In this scenario, just randomize the direction and move on.
      if( the2DSafety >= fEpsilonForWalkOnSpheres ){
	G4ThreeVector returnDir = G4RandomDirection();
	returnDir.setZ(0);
	
	//Debugging
	if( verboseLevel > 5 ){      
	  G4cout << "PSDI Function Point BA | after the FindDirectionToNearbyBoundary returned a funky edge case, the2DSafety is now larger than fEpsilonForWalkOnSpheres." << G4endl;
	  G4cout << "Just randomizing the next direction to " << returnDir.unit() << "." << G4endl;      
	}
      }
      //Otherwise, launch in the direction identified by the safety
      else{	
	returnDir = safetyDir;

	//We note that in the *next* step, running a naive Check2DSafety will yield a distance smaller than that found with this guy. So
	//we need to make a flag that tells the next step about the need to run a special 2D safety based on this edge case being true.
	fNeedSweptSafetyInGetMFP = true;
	

	//Debugging
	if( verboseLevel > 5 ){      
	  G4cout << "PSDI Function Point BB | after the FindDirectionToNearbyBoundary returned a funky edge case, the2DSafety is still smaller than fEpsilonForWalkOnSpheres." << G4endl;
	  G4cout << "Returning a new direction of " << returnDir << " in the global frame." << G4endl;
	}

      }      
    }
    fParticleChange.ProposeMomentumDirection(returnDir.unit());
  }
  //3. If we're not within epsilon, then just purely randomize.
  else{

    //Randomize the final state momentum
    G4ThreeVector returnDir = G4RandomDirection();
    returnDir.setZ(0);
    fParticleChange.ProposeMomentumDirection(returnDir.unit());

    //Debugging
    if( verboseLevel > 5 ){      
      G4cout << "PSDI Function Point BB | safety is larger than epsilon, and we're just randomizing the next direction to " << returnDir.unit() << "." << G4endl;      
    }
  }

  ClearNumberOfInteractionLengthLeft();		// All processes should do this! 
  return &fParticleChange;

}

//Using the information about the track and the safety computed at the track point, we can find a direction to the
//nearby boundary.
G4ThreeVector G4CMPBogoliubovQPRandomWalkTransport::FindDirectionToNearbyBoundary(const G4Track& track, const G4double the2DSafety, G4bool & needToRepeatCalculation, G4bool useSweepForDaughterSafety){

  //Debugging
  if( verboseLevel > 5 ){
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
						     false,
						     useSweepForDaughterSafety);
  //Debugging
  if( verboseLevel > 5 ){
    G4cout << "FDTNB Function Point A | pos: " << track.GetPosition() << ", shiftedPoint: " << shiftedPoint << ", original safety: " << the2DSafety << ", shiftedPoint2Dsafety: " << shiftedPoint2DSafety << G4endl;
  }

  //We now have two points: track.GetPosition() and shiftedPoint, and two safeties: the2DSafety and shiftedPoint2DSafety
  //We can now use these points to identify an angle between the current momDir and the surface normal. For now we'll assume
  //that this angle is in XY, but later (REL) we should come back and fix this to be more plane-agnostic.
  G4double deltaDistToSurface = shiftedPoint2DSafety - the2DSafety;

  //Sometimes our preliminary get2Dsafety functions don't do perfect calculation of the safeties due to the finite
  //granularity of the search in phi. This can occassionally cause the deltaDistToSurface to be slightly larger than the deltaPath,
  //which in principle should never happen. If this does happen, what this implies is that our momentum vector is ALREADY
  //aimed basically directly at the surface. In this case, just return the momentum direction (with sign dependent on the sign of the deltaDistToSurface.
  if( fabs(deltaDistToSurface) > fabs(deltaPath) ){
    if( deltaDistToSurface > 0 ) return momDir;
    else{ return -1*momDir; }
  }

  
  
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
  if( verboseLevel > 5 ){
    G4cout << "FDTNB Function Point B | deltaDistToSurface: " << deltaDistToSurface << G4endl;
    G4cout << "FDTNB Function Point B | deltaPath: " << deltaPath << G4endl;
    G4cout << "FDTNB Function Point B | theta for rotation: " << theta << G4endl;
  }

  //Now, we have the theta, but this theta could be in either direction relative to the surface. So what we'll do is try rotating our momDir vector
  //by theta in both directions, stepping the point forward by each, and re-finding safeties. This calculation is specific to the XY plane and
  //should be generalized (REL).
  G4ThreeVector option1 = momDir.rotateZ(theta);
  G4ThreeVector option2 = momDir.rotateZ(-2*theta);
  momDir.rotateZ(theta);

  //Debugging
  if( verboseLevel > 5 ){
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
  if( verboseLevel > 5 ){
    G4cout << "FDTNB Function Point D | new probe position option 1: " << newPosOption1 << G4endl;
    G4cout << "FDTNB Function Point D | new probe position option 2: " << newPosOption2 << G4endl;
  }

  //Now compute the safeties of our two probe vectors' endpoints to see which is closer.
  G4double option1Safety = G4CMP::Get2DSafety(track.GetStep()->GetPreStepPoint()->GetTouchable(),
					      newPosOption1,
					      track.GetMomentumDirection(),
					      false,
					      useSweepForDaughterSafety);

  //Debugging
  if( verboseLevel > 5 ){
    G4cout << "FDTNB Function Point DA | option 1 safety: " << option1Safety << G4endl;
  }

  
  G4double option2Safety = G4CMP::Get2DSafety(track.GetStep()->GetPreStepPoint()->GetTouchable(),
					      newPosOption2,
					      track.GetMomentumDirection(),
					      false,
					      useSweepForDaughterSafety);
  //Debugging
  if( verboseLevel > 5 ){
    G4cout << "FDTNB Function Point DB | option 2 safety: " << option2Safety << G4endl;
    G4cout << "FDTNB Function Point E | new probe position option 1 safety is " << option1Safety << G4endl;
    G4cout << "FDTNB Function Point E | new probe position option 1 safety x 1e9 is " << option1Safety*1.0e9 << G4endl;
    G4cout << "FDTNB Function Point E | new probe position option 2 safety is: " << option2Safety << G4endl;
    G4cout << "FDTNB Function Point E | new probe position option 2 safety x 1e9 is " << option2Safety*1.0e9 << G4endl;
  }


  //Need to handle an edge case, which happens when the two options' safeties are small enough to be below the picometer G4 tolerance bound
  //and end up being zero. If this is true, bring down the smallerSafety value by a factor and recalculate. This is computationally
  //inefficient but should be rare enough that it hopefully shouldn't matter too much.
  G4double originalOption1Safety = option1Safety;
  G4double originalOption2Safety = option2Safety;
  while( !(option2Safety > 0 && option1Safety > 0) ){
    //    G4ExceptionDescription msg;
    //msg << "Hitting options are both zero edge case" << G4endl;
    //G4Exception("G4CMPBogoliubovQPRandomWalkTransport::FindDirectionToNearbyBoundary", "BogoliubovQPRandomWalkTransport00X",FatalException, msg);
    
    smallerSafety *= 0.5; 
    newPosOption1 = track.GetPosition() + smallerSafety*option1;
    newPosOption2 = track.GetPosition() + smallerSafety*option2;
    option1Safety = G4CMP::Get2DSafety(track.GetStep()->GetPreStepPoint()->GetTouchable(),
				       newPosOption1,
				       track.GetMomentumDirection(),
				       false,
				       useSweepForDaughterSafety);
    option2Safety = G4CMP::Get2DSafety(track.GetStep()->GetPreStepPoint()->GetTouchable(),
				       newPosOption2,
				       track.GetMomentumDirection(),
				       false,
				       useSweepForDaughterSafety);
  }

  
  //If option 1 safety is lower, it means we return option 1 as the direction to the boundary
  G4ThreeVector outputDir;
  if( option1Safety < option2Safety ){
    //if( verboseLevel > 5 ){
      G4cout << "FDNB Function Point F | Direction to the boundary is option 1: " << option1 << G4endl;
      //}
    outputDir = option1;
  }
  else if( option2Safety < option1Safety ){
    //if( verboseLevel > 5 ){
      G4cout << "FDNB Function Point G | Direction to the boundary is option 2: " << option2 << G4endl;
      //}
    outputDir = option2;
  }  
  else{

    //In the edge case where they're both zero, just 
    if( option2Safety == 0 && option1Safety == 0 ){
      G4ExceptionDescription msg;
      msg << "both directions option1 and option2 give the same distance for some reason, and even after a recursive safety recalculation are zero? Should never get here. Option1: " << option1 << ", safety: " << option1Safety << ", Option2: " << option2 << ", safety: " << option2Safety << ", TrackPoint: " << track.GetPosition() << ", the2DSafety (input): " << the2DSafety << ". This seems to be an edge case.";
      G4Exception("G4CMPBogoliubovQPRandomWalkTransport::FindDirectionToNearbyBoundary", "BogoliubovQPRandomWalkTransport00X",FatalException, msg);
    }
    else{
      G4ExceptionDescription msg;
      msg << "both directions option1 and option2 give the same distance for some reason and are NOT zero? Option1: " << option1 << ", safety: " << option1Safety << ", Option2: " << option2 << ", safety: " << option2Safety << ", TrackPoint: " << track.GetPosition() << ", the2DSafety (input): " << the2DSafety << ". This seems to be an edge case.";
      G4Exception("G4CMPBogoliubovQPRandomWalkTransport::FindDirectionToNearbyBoundary", "BogoliubovQPRandomWalkTransport00X",JustWarning, msg);      
    }   
  }

  //Do a penultimate check: move from the current point's location along the output dir by half of the safety and recompute. If the
  //resulting safety is NOT reasonably about 50% of the previous safety, then we need to redo this entire function with the
  //swept daughter safety applied.
  G4double fractionalSafetyDifferenceThreshold = 0.1;
  G4ThreeVector checkPoint = track.GetPosition() + outputDir*the2DSafety*0.5;
  G4double checkedSafety = G4CMP::Get2DSafety(track.GetStep()->GetPreStepPoint()->GetTouchable(),
					      checkPoint,
					      track.GetMomentumDirection(), //I think momDir isn't actually used? Need to consider removing REL
					      false,
					      useSweepForDaughterSafety);
  if( fabs((checkedSafety/the2DSafety) - 0.5) > fractionalSafetyDifferenceThreshold ){
    G4ExceptionDescription msg;
    msg << "When trying to find the direction to the boundary, something seems fishy -- moving in the purported direction of the boundary by half of the safety does not give a new safety that is half of the distance to the boundary. Safety will need to be recomputed with the sweep technique. Pos: " << track.GetPosition() << ", alleged direction to boundary: " << outputDir << ", original safety: " << the2DSafety << ", the checked safety: " << checkedSafety << ", volume at thisPos: " << G4CMP::GetVolumeAtPoint(track.GetPosition()) << " (" << G4CMP::GetVolumeAtPoint(track.GetPosition())->GetName() << G4endl;
    G4Exception("G4CMPBogoliubovQPRandomWalkTransport::FindDirectionToNearbyBoundary", "BogoliubovQPRandomWalkTransport00X",JustWarning, msg);
    needToRepeatCalculation = true;
  }

  //If we pass this and don't need to repeat based on this, do a final check to make sure we're actually crossing a boundary from
  //one volume into another. Sometimes the as-computed direction is to a "phantom boundary" because G4 doesn't calculate all "quick"
  //safeties as accurately as possible.
  if( needToRepeatCalculation == false ){
    G4ThreeVector thisPos = track.GetPosition();
    G4ThreeVector posOstensiblyOverBoundary(0,0,0);

    //Edge case: we're close enough that the two options are originally effectively zero safety. Here let's just artificially inflate our nudge by a bit
    if( originalOption1Safety == 0 && originalOption2Safety == 0 ){
      posOstensiblyOverBoundary = track.GetPosition() + the2DSafety*outputDir*fBoundaryFudgeFactor*2;
    }
    //Most common
    else{
      posOstensiblyOverBoundary = track.GetPosition() + the2DSafety*outputDir*fBoundaryFudgeFactor;
    }


    if( G4CMP::GetVolumeAtPoint(thisPos) == G4CMP::GetVolumeAtPoint(posOstensiblyOverBoundary) ){
      G4ExceptionDescription msg;
      msg << "When trying to find the direction to the boundary, we seem to be calculating a direction vector that satisfies our first, fractionalSafetyDifference criterion, but not our phantom boundary condition. (This may sometimes happen near corners, where our direction-finding algorithm's math breaks down. At position: " << track.GetPosition() << ", alleged direction to boundary: " << outputDir << ", original safety: " << the2DSafety << ", volume at thisPos: " << G4CMP::GetVolumeAtPoint(thisPos) << " (" << G4CMP::GetVolumeAtPoint(thisPos)->GetName() << ", volume at the position ostensibly over the boundary: " << G4CMP::GetVolumeAtPoint(posOstensiblyOverBoundary) << G4endl;
      G4Exception("G4CMPBogoliubovQPRandomWalkTransport::FindDirectionToNearbyBoundary", "BogoliubovQPRandomWalkTransport00X",JustWarning, msg);
      needToRepeatCalculation = true;
    }
  }
  
  //One ACTUALLY final check: compute the volume at the current point and the volume at a point that is one checkedSafety away

  
  return outputDir;
}
					     



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4CMPBogoliubovQPRandomWalkTransport::GetContinuousStepLimit(
                                       const G4Track& track,
                                       G4double previousStepSize,
                                       G4double currentMinimalStep,
                                       G4double& currentSafety)
{
  //Debugging
  if( verboseLevel > 5 ){
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
  if( verboseLevel > 5 ){
    G4cout << "---------- G4CMPBogoliubovQPRandomWalkTransport::GetMeanFreePath ----------" << G4endl;
    G4cout << "GMFP Function Point A | track volume: " << track.GetVolume()->GetName() << G4endl;
    G4cout << "GMFP Function Point A | previousStepSize: " << previousStepSize << G4endl;
    G4cout << "GMFP Function Point A | momentum direction: " << track.GetMomentumDirection() << G4endl;
    G4cout << "GMFP Function Point A | position: " << track.GetPosition() << G4endl;
    G4cout << "GMFP Function Point A | global time: " << track.GetGlobalTime() << G4endl;
    G4cout << "GMFP Function Point A | forcing condition: " << *condition << ", compared to Forced: " << Forced << " and NotForced: " << NotForced << G4endl;
  }

  //I think this is needed here -- we're overriding the discrete GPIL functions in the base G4VContinuousDiscrete class,
  //but since we haven't set the condition variable anywhere else, I think this ends up taking on a garbage value unless
  //we set it. Here we set it to not forced, which makes this class's discrete bit compete on equal footing with the other
  //discrete processes (as it should)
  *condition = NotForced;
  
  //Reset some of the quantities that should be reset each step
  fOutgoingSurfaceTangent1 = G4ThreeVector(0,0,0);
  fOutgoingSurfaceTangent2 = G4ThreeVector(0,0,0);
  fQPIsStuck = false;
  
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
  if( verboseLevel > 5 ){
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
    //pushes us across a corner back into World, we'll have an issue. Hopefully this only happens in very "badly-behaved" geometries (like 20-
    //pointed stars or something with highly acute internal corners), since the momentum should in principle be the normal vector based on
    //the boundary physics.
    
    //Debugging
    if( verboseLevel > 5 ){      
      G4cout << "GMFP Function Point C | In a turnaround step. Killing the transport GPIL." << G4endl;
    }    
    fTrackOnBoundary = true;
    isActive = false;
    return DBL_MAX;
  }

  //We can be on a boundary and still be "active" if we're in the volume into which we're moving.
  if( theStatus == fGeomBoundary ){
    
    //Debugging
    if( verboseLevel > 5 ){      
      G4cout << "GMFP Function Point D | On a boundary but not in a turnaround step." << G4endl;
    }
    fTrackOnBoundary = true;
  }

  //Debugging
  if( verboseLevel > 5 ){
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

      //Here we need to split into "normal" and "QP-is-stuck" scenarios. Have to do this here, and not later, so that we can
      //properly enable the GPIL race with other processes in either case. First, establish the surface norm. From the boundary
      //process, this should just be the momentum direction at the pre-step point (track point). Then check to see if the QP is stuck
      G4ThreeVector surfaceNorm = track.GetMomentumDirection();      
      UpdateBoundaryHistory(track.GetTrackID(),track.GetPosition(),surfaceNorm);
      std::tuple<G4bool,G4ThreeVector,G4ThreeVector,G4ThreeVector,G4ThreeVector> qpIsStuck_norm1_norm2_pos1_pos2 = CheckForStuckQPs();
      G4bool qpIsStuck = std::get<0>(qpIsStuck_norm1_norm2_pos1_pos2);
      

      //If our QP is not stuck, then we get an easy win. Compute the 2D safety normally
      if( !qpIsStuck ){

	//Debugging
	if( verboseLevel > 5 ){
	  G4cout << "GMFP Function Point EB | QP Is not stuck." << G4endl;
	}

	//If the post-step do it of last step suggests that we need a swept safety to compensate for G4's incorrectness, then run a swept
	//2D safety from the boundary. I think this usually shouldn't run.
	if(!fNeedSweptSafetyInGetMFP){
	  the2DSafety = G4CMP::Get2DSafety(track.GetStep()->GetPreStepPoint()->GetTouchable(),
					   track.GetPosition(),
					   track.GetMomentumDirection(),
					   true,
					   false,
					   surfaceNorm );
	}
	else{	 
	  std::pair<G4double,G4ThreeVector> the2DSafetyAndDir = G4CMP::Get2DSafetyWithDirection(track.GetStep()->GetPreStepPoint()->GetTouchable(),
												track.GetPosition(),
												track.GetMomentumDirection(),
												true,
												surfaceNorm);
	  the2DSafety = the2DSafetyAndDir.first;

	  G4ExceptionDescription msg;
	  msg << "In GetMFP We're somehow on a boundary and also have triggered the fNeedSweptSafetyInGetMFP. What is happening? (In principle, this isn't a deathknell-- I'm just killing the code here to see where this even is triggered.";
	  G4Exception("G4CMPBogoliubovQPRandomWalkTransport::GetMeanFreePath", "BogoliubovQPRandomWalkTransport002",FatalException, msg);
	}
      }
      //If the QP IS stuck, then we have some work to do. Use the norms and positions returned from CheckForStuckQPs to compute
      //an angular range away from the corner over which we can compute a safety.
      else{

	G4ThreeVector norm1 = std::get<1>(qpIsStuck_norm1_norm2_pos1_pos2);
	G4ThreeVector pos1 = std::get<3>(qpIsStuck_norm1_norm2_pos1_pos2);
	G4ThreeVector norm2 = std::get<2>(qpIsStuck_norm1_norm2_pos1_pos2);
	G4ThreeVector pos2 = std::get<4>(qpIsStuck_norm1_norm2_pos1_pos2);

	//From these normal vectors, we should be able to find a set of directions over which to do a constrained 2D safety and 
	//subsequent ejection of the QP out from the corner
	G4ThreeVector cornerLocation(0,0,0);
	std::pair<G4ThreeVector,G4ThreeVector> outgoingSurfaceTangents = FindSurfaceTangentsForStuckQPEjection(norm1,pos1,norm2,pos2,cornerLocation);
	fOutgoingSurfaceTangent1 = outgoingSurfaceTangents.first;
	fOutgoingSurfaceTangent2 = outgoingSurfaceTangents.second;
	fQPIsStuck = true;

	//Debugging
	if( verboseLevel > 5 ){
	  G4cout << "GMFP Function Point EB | QP Is stuck!" << G4endl;
	  G4cout << "GMFP Function Point EB | norm 1: " << norm1 << G4endl;
	  G4cout << "GMFP Function Point EB | norm 2: " << norm2 << G4endl;
	  G4cout << "GMFP Function Point EB | pos 1: " << pos1 << G4endl;
	  G4cout << "GMFP Function Point EB | pos 2: " << pos2 << G4endl;
	  G4cout << "GMFP Function Point EB | fOutgoingSurfaceTangent1: " << fOutgoingSurfaceTangent1 << G4endl;
	  G4cout << "GMFP Function Point EB | fOutgoingSurfaceTangent2: " << fOutgoingSurfaceTangent2 << G4endl;
	  G4cout << "GMFP Function Point EB | cornerLocation: " << cornerLocation << G4endl;
	}
	
	
	//Now use the outgoing surface tangents to compute a constrained 2D safety
	G4double constrained2DSafety = G4CMP::Get2DSafety(track.GetStep()->GetPreStepPoint()->GetTouchable(),
							  track.GetPosition(),
							  track.GetMomentumDirection(),
							  true,
							  false,
							  G4ThreeVector(0,0,0), //The norm should NOT be well-defined if we're stuck in a corner
							  outgoingSurfaceTangents.first,
							  outgoingSurfaceTangents.second);	
	the2DSafety = constrained2DSafety;

	//Debugging
	if( verboseLevel > 5 ){
	  G4cout << "GMFP Function Point EC | Constrained 2D safety: " << constrained2DSafety << G4endl;
	}
      }
    }
    //Bulk case: simpler
    else{

      //If we don't need a swept safety as determined by the last step, then run the normal (faster) one
      if(!fNeedSweptSafetyInGetMFP){      
	the2DSafety = G4CMP::Get2DSafety(track.GetStep()->GetPreStepPoint()->GetTouchable(),
					 track.GetPosition(),
					 track.GetMomentumDirection(),
					 false);
      }
      //If we need the swept safety to compensate for G4 handling a daughter safety wrong, then run the swept safety.
      else{
	std::pair<G4double,G4ThreeVector> the2DSafetyAndDir = G4CMP::Get2DSafetyWithDirection(track.GetStep()->GetPreStepPoint()->GetTouchable(),
											      track.GetPosition(),
											      track.GetMomentumDirection(),
											      true);
	the2DSafety = the2DSafetyAndDir.first;
	G4cout << "the2DSafety after doing a swept thing in GetMFP: " << the2DSafety << G4endl;
      }
    }
    f2DSafety = the2DSafety;
    
    //Calculate the energy dependent diffusion constant
    G4double E_ratio = fGapEnergy/energy;
    G4double E_ratio2 = pow(E_ratio,2.0);
    fDiffConst = fDn*sqrt(1-E_ratio2);

    //Debugging
    if( verboseLevel > 5 ){      
      G4cout << "GMFP Function Point F | Diffusion constant (energy-adjusted): " << fDiffConst << G4endl;
      G4cout << "GMFP Function Point F | the2DSafety: " << the2DSafety << G4endl;
    }
      
    
    //Using this diffusion constant, draw a time to reach the boundary
    //established by the safety. For now we use something simple and wrong to test,
    //taking the 2D safety to be the sigma of the diffusion broadening
    //but need the true expression from Grebenkov
    G4double timeStepToBoundary = SampleTimeStepFromFirstPassageDistribution(the2DSafety);    
    fTimeStepToBoundary = timeStepToBoundary;

    //Debugging
    if( verboseLevel > 5 ){      
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
    //REL 4/1/2025 Now that we have TOFP in, I'm not sure if this matters any more... Should check this. Actually, think that the
    //important thing here is actually the relationship between diffusion constant and velocity, which is I think
    //somewhat shaky since we keep a constant velocity but allow for energy-dependent diffusion. So actually yeah I think
    //we have to keep this for generality.
    fVerySmallStep = false;
    if( thisMFP < the2DSafety ){
      fVerySmallStep = true;
      thisMFP = the2DSafety*fBoundaryFudgeFactor; //Fudge factor needs to go here as well because AlongStep GPIL uses PhysicalStep, which comes from PostStepGPIL (this is so that we can enable Transportation to trigger/win the AlongStep race)
      
      //G4ExceptionDescription msg;
      //msg << "We're triggering fVerySmallStep even though we're using the time-of-first-passage formalism.";
      //G4Exception("G4CMPBogoliubovQPRandomWalkTransport::GetMeanFreePath", "BogoliubovQPRandomWalkTransport002",FatalException, msg);
    }
    
    //Now, as a last measure, we set the number of interaction lengths left for this process to 1
    //and return a mean free path equal to this MFP, so that the only way that another process will
    //win over this in the discrete GPIL race is if the time to that process is less than the time
    //to this one.
    theNumberOfInteractionLengthLeft = 1;

    //Debugging
    if( verboseLevel > 5 ){      
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
  if( verboseLevel > 5 ){
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
    if( verboseLevel > 5 ){
      G4cout << "UMFPFLC Function Point B | Step length associated with this is " << aTrack.GetStep()->GetStepLength() << G4endl;
      G4cout << "UMFPFLC Function Point B | Successfully changed over to a new lattice for process " << this->GetProcessName() << G4endl;
    }
        
    //REL noting that if physical lattices are not 1:1 with volumes, something may get broken here... Should check a scenario of segmented SC...    
    this->LoadDataForTrack(&aTrack);
    return true;
    
  }

  //Debugging
  if( verboseLevel > 5 ){
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
  if( verboseLevel > 5 ){
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
    return 0; //Shouldn't ever get here.
  }
}

//This samples a dimensionless time step using A/R technique. 
G4double G4CMPBogoliubovQPRandomWalkTransport::SampleDimensionlessTimeStepUsingAcceptanceRejectionTechnique()
{
  if( verboseLevel > 5 ){
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
    if( verboseLevel > 5 ){
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

//This looks at the storage objects and uses them to understand if our track hasn't really moved
//much in the last several boundary scatters. The outputs are whether a QP is stuck and then two
//vectors that give a span over which to search for a new 2D safety.
std::tuple<G4bool,G4ThreeVector,G4ThreeVector,G4ThreeVector,G4ThreeVector> G4CMPBogoliubovQPRandomWalkTransport::CheckForStuckQPs(){
  
  //For now, just check for stuck QPs in a corner geometry
  return CheckForStuckQPsInCorner();
}

//For a specific geometry where we are explicitly dealing with a corner (doesn't have to be 90 degrees, but it has to
//be discontinuous). General strategy:
std::tuple<G4bool,G4ThreeVector,G4ThreeVector,G4ThreeVector,G4ThreeVector> G4CMPBogoliubovQPRandomWalkTransport::CheckForStuckQPsInCorner(){

  //Debugging
  if( verboseLevel > 5 ){
    G4cout << "---------- G4CMPBogoliubovQPRandomWalkTransport::CheckForStuckQPsInCorner() ----------" << G4endl;
    G4cout << "CFSQIC Function Point A | Starting to check for stuck QPs in a corner." << G4endl;
  }
  
  
  //output variables
  G4bool qpIsStuck = false;
  G4ThreeVector outNorm0(0,0,0);
  G4ThreeVector outNorm1(0,0,0);
  G4ThreeVector outPos0(0,0,0);
  G4ThreeVector outPos1(0,0,0);

  //Some storage variables
  G4double avgx = 0;
  G4double avgx2 = 0;
  G4double avgy = 0;
  G4double avgy2 = 0;
  std::vector<G4ThreeVector> uniqueNorms;

  //First, if the length of the boundary history is not the max length, it means we haven't been running for long
  //enough to be stuck. Return false
  if( fBoundaryHistory.size() < fMaxBoundaryHistoryEntries ){

    //Debugging
    if( verboseLevel > 5 ){
      G4cout << "CFSQIC Function Point AB | fBoundaryHistory.size() is less than fMaxBoundaryHistoryEntries. Aborting stuck QP check." << G4endl;
    }    
    std::tuple<G4bool,G4ThreeVector,G4ThreeVector,G4ThreeVector,G4ThreeVector> output(qpIsStuck,outNorm0,outNorm1,outPos0,outPos1);
    return output;
  }
  
  //Find the most recent boundary history norm and position, and which volume is in the "forward" direction
  G4double epsilonDisplacement = 1 * CLHEP::nm;
  G4ThreeVector currentPosition = fBoundaryHistory[fBoundaryHistory.size()-1].first;
  G4ThreeVector currentNorm = fBoundaryHistory[fBoundaryHistory.size()-1].second;
  G4ThreeVector displacedPosition = currentPosition + currentNorm*epsilonDisplacement;
  G4VPhysicalVolume * volumeAtPoint = G4CMP::GetVolumeAtPoint(displacedPosition);

  //Debugging
  if( verboseLevel > 5 ){
    G4cout << "CFSQIC Function Point AC | current position: " << currentPosition << G4endl;
    G4cout << "CFSQIC Function Point AC | current norm: " << currentNorm << G4endl;
    G4cout << "CFSQIC Function Point AC | volumeAtPoint: " << volumeAtPoint->GetName() << G4endl;
  }


  
  //Loop over the other boundary history info (here we include the "current" one but a clause that will continue if it
  //sees it.
  std::vector<G4ThreeVector> goodOtherNorms;
  std::vector<G4ThreeVector> goodOtherPositions;
  for( int iB = 0; iB < fBoundaryHistory.size(); ++iB ){

    //First, compute the standard deviation of the positions. This should be done regardless of the
    //norm, since it will give us a sense of how tightly clustered the boundary hits are.
    double x = fBoundaryHistory[iB].first.getX();
    double y = fBoundaryHistory[iB].first.getY();
    avgx += x;
    avgy += y;
    avgx2 += x*x;
    avgy2 += y*y;
    
    //Now Determine if the currentNorm is collinear with this historical norm. If it is, continue: don't want
    //things that are parallel with the currentNorm
    G4ThreeVector thisNorm = fBoundaryHistory[iB].second;
    if( fabs(thisNorm.dot(currentNorm)) > fDotProductDefiningUniqueNorms ){
      continue;
    }

    //Since we now have a norm that is not collinear, compute +/- eps distances from this position along this norm,
    //and determine if it's in the same volume as the current norm's next volume. This should leave us with
    //a set of norms and points corresponding to the other positions that point into the same volume as the "main" one
    G4ThreeVector plusEps = fBoundaryHistory[iB].first + epsilonDisplacement * fBoundaryHistory[iB].second;
    G4ThreeVector minusEps = fBoundaryHistory[iB].first - epsilonDisplacement * fBoundaryHistory[iB].second;
    G4VPhysicalVolume * volPlusEps = G4CMP::GetVolumeAtPoint(plusEps);
    G4VPhysicalVolume * volMinusEps = G4CMP::GetVolumeAtPoint(minusEps);
    if( volPlusEps == volumeAtPoint ){
      goodOtherNorms.push_back(fBoundaryHistory[iB].second);
      goodOtherPositions.push_back(fBoundaryHistory[iB].first);
      if( verboseLevel > 5 ){
	G4cout << "CFSQIC Function Point AD | PlusEps case. Goodothernorm: " << fBoundaryHistory[iB].second << ", pos: " << fBoundaryHistory[iB].first << G4endl;
      }
    }
    if( volMinusEps == volumeAtPoint ){
      goodOtherNorms.push_back(-1*fBoundaryHistory[iB].second);
      goodOtherPositions.push_back(fBoundaryHistory[iB].first);
      if( verboseLevel > 5 ){
	G4cout << "CFSQIC Function Point AE | MinusEps case. Goodothernorm: " << fBoundaryHistory[iB].second << ", pos: " << fBoundaryHistory[iB].first << G4endl;
      }
    }
  }
  avgx /= fBoundaryHistory.size();
  avgy /= fBoundaryHistory.size();
  avgx2 /= fBoundaryHistory.size();
  avgy2 /= fBoundaryHistory.size();
  G4double sigmax = pow(avgx2 - avgx*avgx,0.5);
  G4double sigmay = pow(avgy2 - avgy*avgy,0.5);

  //Debugging
  if( verboseLevel > 5 ){
    G4cout << "CFSQIC Function Point B | avgx: " << avgx << ", avgx2: " << avgx2 << G4endl;
    G4cout << "CFSQIC Function Point B | sigmaX: " << sigmax << ", sigmaY: " << sigmay << G4endl;
    G4cout << "CFSQIC Function Point B | Length of goodOtherNorms: " << goodOtherNorms.size() << G4endl;
  }
  
  //Next, check to see that the sigmaX and sigmaY are both below a threshold. If they're not, then return that we're not stuck.
  if( !(sigmax < fStuckInCornerThreshold && sigmay < fStuckInCornerThreshold ) ){
    std::tuple<G4bool,G4ThreeVector,G4ThreeVector,G4ThreeVector,G4ThreeVector> output(qpIsStuck,outNorm0,outNorm1,outPos0,outPos1);

    //Debugging
    if( verboseLevel > 5 ){
      G4cout << "CFSQIC Function Point BA | Looks like we're not stuck. Either sigmax or sigmay is large enough to be not stuck." << G4endl;
    }
    return output;    
  }

  //Next, check to see if we have any remaining good other norms. If, for example, we're on a curved surface, that number may be zero. We'll
  //start by just saying we're not stuck in that case.
  if( goodOtherNorms.size() == 0 ){
    std::tuple<G4bool,G4ThreeVector,G4ThreeVector,G4ThreeVector,G4ThreeVector> output(qpIsStuck,outNorm0,outNorm1,outPos0,outPos1);

    //Debugging
    if( verboseLevel > 5 ){
      G4cout << "CFSQIC Function Point BB | Looks like we have a cluster of close points but have zero good other norms, which may occur if we're on a curved surface. For now we'll say we're not stuck." << G4endl;
    }    
    return output;
  }
  
  outPos0 = currentPosition;
  outNorm0 = currentNorm;

  //Debugging
  if( verboseLevel > 5 ){
    for( int iN = 0; iN < goodOtherNorms.size(); ++iN ){
      G4cout << "CFSQIC Function Point BC | Good other norm: " << goodOtherNorms[iN] << ", pos: " << goodOtherPositions[iN] << G4endl;
    }
  }
    
  
  //So we're boxed into a corner. Now check the good positions list to see if there are more than one unique vector that points to this next
  //volume. In the case of curved surfaces, this may very well be true.
  std::vector<G4ThreeVector> uniqueGoodNorms;
  std::vector<G4ThreeVector> uniqueGoodPositions;
  uniqueGoodNorms.push_back(goodOtherNorms[0]);
  uniqueGoodPositions.push_back(goodOtherPositions[0]);

  //Loop over all of the other good norms and remove if not unique
  for( int iN = 0; iN < goodOtherNorms.size(); ++iN ){
    G4bool notUnique = false;
    for( int iG = 0; iG < uniqueGoodNorms.size(); ++iG ){
      if( fabs(goodOtherNorms[iN].dot(uniqueGoodNorms[iG])) > fDotProductDefiningUniqueNorms ){
	notUnique = true;
	break;
      }
    }
    if( !notUnique ){
      uniqueGoodNorms.push_back(goodOtherNorms[iN]);
      uniqueGoodPositions.push_back(goodOtherPositions[iN]);
    }
  }
  
  //Debugging
  if( verboseLevel > 5 ){
    G4cout << "CFSQIC Function Point C | Unique Good *current* Norm: " << currentNorm << " at position: " << currentPosition << G4endl;
    G4cout << "CFSQIC Function Point C | Number of unique good norms: " << uniqueGoodNorms.size() << ", number of corresp. positions: " << uniqueGoodPositions.size() << G4endl;
    for( int iU = 0; iU < uniqueGoodNorms.size(); ++iU ){
      G4cout << "CFSQIC Function Point C | Unique Good *other* Norm: " << uniqueGoodNorms[iU] << " at position: " << uniqueGoodPositions[iU] << G4endl;
    }
  }
    
  
  //Now we should have a set of unique good norms, pointing into the volume of interest. For now, if we're really in a true
  //discontinuous corner, the length of this unique good norms vector should never be larger than one. However, in some
  //more general geometries it may -- we'll have to figure out what we want to do about that later. (REL PLACE FOR MORE GENERALITY)
  if( uniqueGoodNorms.size() != 1 ){
    G4ExceptionDescription msg;
    msg << "The number of unique good norms at this corner is somehow not equal to 2. This is probably because we're looking at a curved geometry. We don't have logic to handle these yet. Apologies.";
    G4Exception("G4CMPBogoliubovQPRandomWalkTransport::CheckForStuckQPsInCorner", "BogoliubovQPRandomWalkTransport00X",FatalException, msg);
  }
  
  //Now package the good norms and positions into output
  outPos1 = uniqueGoodPositions[0];
  outNorm1 = uniqueGoodNorms[0];
  qpIsStuck = true;
  std::tuple<G4bool,G4ThreeVector,G4ThreeVector,G4ThreeVector,G4ThreeVector> output(qpIsStuck,outNorm0,outNorm1,outPos0,outPos1);
  return output;
}



//This takes the track/step information and pushes it back into our storage objects that keep track
//of the last few boundary steps
void G4CMPBogoliubovQPRandomWalkTransport::UpdateBoundaryHistory(G4int trackID, G4ThreeVector preStepPos, G4ThreeVector preStepNorm){

  //Debugging
  if( verboseLevel > 5 ){
    G4cout << "---------- G4CMPBogoliubovQPRandomWalkTransport::UpdateBoundaryHistory() ----------" << G4endl;
  }
  
  //Check the track id to make sure we're looking at the right track. If they don't match,
  //reset our boundary history and push our vector back
  if( trackID != fBoundaryHistoryTrackID ){
    fBoundaryHistory.clear();
    fBoundaryHistoryTrackID = trackID;
    std::pair<G4ThreeVector,G4ThreeVector> posNormPair(preStepPos,preStepNorm);    
    fBoundaryHistory.push_back(posNormPair);
  }
  //If we're on the same particle, then push back our boundary history vectors until they
  //hit a max length, and then start kicking old boundary interactions out, FIFO  
  else{
    std::pair<G4ThreeVector,G4ThreeVector> posNormPair(preStepPos,preStepNorm);    
    if( fBoundaryHistory.size() < fMaxBoundaryHistoryEntries ){
      fBoundaryHistory.push_back(posNormPair);      
    }
    else{
      fBoundaryHistory.erase(fBoundaryHistory.begin());
      fBoundaryHistory.push_back(posNormPair);
    }    
  }

  //Debugging
  if( verboseLevel > 5 ){
    for( int iB = 0; iB < fBoundaryHistory.size(); ++iB ){
      G4cout << "UBH Function Point A | Boundary position: " << fBoundaryHistory[iB].first.getX() << ", " << fBoundaryHistory[iB].first.getY() << G4endl;
    }
  }

  
}

//This takes the positions and norms from the "am I stuck" function and uses them to find the outgoing tangent vectors from a corner
//Note that this function is not yet agnostic to the plane: it's only in XY for now. This can be plausibly made more general.
std::pair<G4ThreeVector,G4ThreeVector> G4CMPBogoliubovQPRandomWalkTransport::FindSurfaceTangentsForStuckQPEjection(G4ThreeVector norm1,
														   G4ThreeVector pos1,
														   G4ThreeVector norm2,
														   G4ThreeVector pos2,
														   G4ThreeVector & cornerLocation)
{
  //Debugging
  if( verboseLevel > 5 ){
    G4cout << "---------- G4CMPBogoliubovQPRandomWalkTransport::FindSurfaceTangentsForStuckQPEjection() ----------" << G4endl;
    G4cout << "FSTFSQE Function Point A | Starting to find surface tangents." << G4endl;
  }

  G4double epsilonCornerCalculation = 0.01*nm; //Should be pretty small -- this all should be exact
  
  //Compute an in-plane vector using the difference between the positions, and then cross with the norm to get the out-of-plane vector  
  G4ThreeVector inPlane(pos1.getX()-pos2.getX(),pos1.getY()-pos2.getY(),pos1.getZ()-pos2.getZ());
  G4ThreeVector outOfPlane = (inPlane.cross(norm1)).unit();

  //Handle floating point errors that might push the norm into slight non-orthogonality. For now this is not geometry-orientation-agnostic,
  //but is a needed sanity check.
  if( outOfPlane.getX() != 0 || outOfPlane.getY() != 0 ){
    if( fabs(outOfPlane.getX()) < 1e-10 ){ outOfPlane.setX(0); }
    if( fabs(outOfPlane.getY()) < 1e-10 ){ outOfPlane.setY(0); }
    if( outOfPlane.getX() != 0 || outOfPlane.getY() != 0 ){
      G4ExceptionDescription msg;
      msg << "Currently, we are only working in XY. After adjustment for reasonable floating point errors, the out-of-plane vector is still somehow not purely along z.";
      G4Exception("G4CMPBogoliubovQPRandomWalkTransport::FindSurfaceTangentsForStuckQPEjection", "BogoliubovQPRandomWalkTransport00X",FatalException, msg);
    }
  }
  
  
  //Now use the outOfPlane vector with the norms to get tangent vectors (with still-to-be-defined directions relative to the corner)
  G4ThreeVector tangVect1 = (outOfPlane.cross(norm1)).unit();
  G4ThreeVector tangVect2 = (outOfPlane.cross(norm2)).unit();

  //Debugging
  if( verboseLevel > 5 ){
    G4cout << "FSTFSQE Function Point B | Input norm1: " << norm1 << G4endl;
    G4cout << "FSTFSQE Function Point B | Input norm2: " << norm2 << G4endl;
    G4cout << "FSTFSQE Function Point B | Additional In-plane vector: " << inPlane << G4endl;
    G4cout << "FSTFSQE Function Point B | Out-of-plane vector: " << outOfPlane << G4endl;
    G4cout << "FSTFSQE Function Point B | Simple TangVect1: " << tangVect1 << G4endl;
    G4cout << "FSTFSQE Function Point B | Simple TangVect2: " << tangVect2 << G4endl;
  }

  //Now from these, compute the location of the corner. This is fully vectorial, and fully general.
  //Here we use a prescription and notation in songho.ca/math/line/line.html where we have two lines, P1 + t * v and Q1 + s * u
  //where we take pos1 = P1, pos2 = Q, tangVect1 = v, tangVect2 = u.
  //We first compute the parameter t at which those two lines are equal to each other:
  G4ThreeVector numPart1 = (pos2-pos1).cross(tangVect2);
  G4ThreeVector numPart2 = tangVect1.cross(tangVect2);
  G4double num = numPart1.dot(numPart2);
  G4ThreeVector denomPart1 = tangVect1.cross(tangVect2);
  G4double denom = denomPart1.dot(denomPart1);
  G4double t = num/denom;

  //Sanity check: can recalculate with s and find

  //Calculating t and s
  //if( verboseLevel > 5 ){
  //  G4cout << "FSTFSQE Function Point B | t_param: " << t_param << G4endl;
  //  G4cout << "FSTFSQE Function Point B | s_param: " << s_param << G4endl;
  //}
  
  //Compute the corner location
  G4ThreeVector corner1 = pos1 + tangVect1*t;
  //if( (corner1-corner2).mag() > epsilonCornerCalculation ){
  //  G4ExceptionDescription msg;
  //  msg << "Somehow our corner calculation has failed. Corner 1: " << corner1 << ", and Corner 2: " << corner2 << ".";
  //  G4Exception("G4CMPBogoliubovQPRandomWalkTransport::FindSurfaceTangentsForStuckQPEjection", "BogoliubovQPRandomWalkTransport00X",FatalException, msg);
  //}
  cornerLocation = corner1;

  //Now compute a vector from the corner to each position
  G4ThreeVector finalTangVect1 = (pos1-corner1).unit();
  G4ThreeVector finalTangVect2 = (pos2-corner1).unit();
  std::pair<G4ThreeVector,G4ThreeVector> output(finalTangVect1,finalTangVect2);
  return output;
}



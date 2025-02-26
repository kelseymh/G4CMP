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
  fBoundaryFudgeFactor = 1.0001;

  
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
//This is the "stretch" version
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
  
  //At this point, the currentMinimalStep is the one that has won the race. We'll take this one 

  
  //Set the path length and pre-diffusion path length to the length suggested by discrete process race winner
  fPathLength = currentMinimalStep;
  fPreDiffusionPathLength = currentMinimalStep;
  G4cout << "---> REL/EY: fPathLength in ASGPIL: " << currentMinimalStep << G4endl;
  G4cout << "---> REL/EY: fPreDiffusionPathLength in ASGPIL: " << currentMinimalStep << G4endl;
  G4cout << "---> REL/EY: velocity in ASGPIL: " << track.GetVelocity() << G4endl;

  
  //Get energy and velocity of track
  G4double energy = track.GetKineticEnergy();
  G4double velocity = track.GetVelocity();
  G4ThreeVector momentumDir = track.GetMomentumDirection();

  //If we're in a turnaround step, then kill
  if( isActive == false ){
    G4cout << "---> REL/EY ASGIPL: In a turnaround step. Killing the transport GPIL." << G4endl;
    return DBL_MAX;
  }
  
  
  /*
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

  */


  
  //3. The diffusion constant and other superconductor-specific information are accessed via the
  //   SCUtils class that this inherits from. No need to actually grab it from somewhere: it's already
  //   a data member.
  G4cout << "---> REL/EY: The gap energy, drawn from SCUtils, is: " << fGapEnergy << G4endl;
  G4cout << "---> REL/EY: The Dn, drawn from SCUtils, is: " << fDn << ", and the fDiffConst: " << fDiffConst << G4endl;
  G4cout << "---> REL/EY: The Teff, drawn from SCUtils, is: " << fTeff << G4endl;  
  G4cout << "---> REL/EY: Gap energy: " << fGapEnergy << ", energy: " << energy << ", DN: " << fDn << G4endl; 

  //4. Given this, compute a step length based on the time step.
  if ((energy>=fGapEnergy) && (fDn)>0){ 
    isActive = true;
    *selection = NotCandidateForSelection;

    //At this point, we just need to return a diffusion-displacement path length. This takes on two forms:
    //1. When the 2D safety in this process limits the step, we need to return exactly this radius.
    //2. When some other process limits the step, we need to return a radius corresponding to that time.
    //   The selected radius for that time MUST be less than the boundary-limited radius, because otherwise the
    //   particle would have passed through that radius and therefore the boundary-limited time to first crossing
    //   would be longer than that computed in option 1. I think we maybe just draw from a gaussian until we
    //   get a radius that's less than the boundary number?

    //Boundary-limited case
    G4cout << "--->REL/EY: currentMinimalStep/velocity: " << currentMinimalStep/velocity << ", fTimeStepToBoundary: " << fTimeStepToBoundary << G4endl;
    G4cout << "--->REL/EY: fTimeStepToBoundary: " << fTimeStepToBoundary << G4endl;

    //Check to see if we're in a boundary-limited case
    double timeTolerance = 1E-10; //For floating point errors    
    if( (fabs(currentMinimalStep/velocity - fTimeStepToBoundary) < timeTolerance) || fVerySmallStep ){
      fPathLength = f2DSafety*fBoundaryFudgeFactor; //Fudge factor so that we can trigger transportation when we are pointed at the boundary exactly
      fTimeStep = fTimeStepToBoundary;
      G4cout << "--->REL/EY: In ASGPIL: Looks like the boundary-limited case applies here. Returning fPathLength = f2DSafety = " << fPathLength << G4endl;
    }
    //Case limited by another physics process
    else{
      fTimeStep = currentMinimalStep / velocity;
      G4double sigma1D = pow(2*fDiffConst*fTimeStep,0.5);
      G4cout << "--->REL/EY: In ASGPIL: Sigma1D: " << sigma1D << " from fDiffConst: " << fDiffConst << G4endl;
      G4RandGauss* gauss_dist = new G4RandGauss(G4Random::getTheEngine(),0.0,sigma1D);
      do{
	double gauss_dist_x = fabs(gauss_dist->fire());
	double gauss_dist_y = fabs(gauss_dist->fire());
	fPathLength = pow(gauss_dist_x*gauss_dist_x + gauss_dist_y*gauss_dist_y,0.5);
      }
      while( fPathLength > f2DSafety );
      
      G4cout << "--->REL/EY: In ASGPIL: Looks like a different discrete process wins the GPIL race. Returning fPathLength = " << fPathLength << G4endl;
    }
    
    G4cout << "--->REL/EY: Successfully returning AlongStepGPIL (diffusion-folded) fPathLength of " << fPathLength << G4endl;
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


/*
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// This is the "working" version
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
*/


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
  return mfp;

  //What was there before
  //*condition = NotForced;
  //return DBL_MAX;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// "Stretch" version
G4VParticleChange* G4CMPBogoliubovQPRandomWalkTransport::AlongStepDoIt(const G4Track& track, const G4Step& step) {

  G4cout << "REL-- InRandomWalkTransport::AlongStepDoIt() A" << G4endl;
  G4cout << "REL: firstlook status1: " << track.GetStep()->GetPreStepPoint()->GetStepStatus() << G4endl;
  G4cout << "REL: firstlook status2: " << step.GetPreStepPoint()->GetStepStatus() << G4endl;
  G4cout << "REL: firstlook status3: " << step.GetPostStepPoint()->GetStepStatus() << G4endl;

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
    //If we hit the wall, we should have set things up so that:
    //1. The discrete process leading to a wall hit is the boundary-limited one set here
    //2. When we hit the wall, we are hitting it mostly "at the right time," so the above code (without the "relocate" bit) still applies
    else if( nextStepLength != kInfinity ){ //We're on a surface
      if( step.GetPostStepPoint()->GetStepStatus() != fGeomBoundary ){
	G4ExceptionDescription msg;
	msg << "Somehow the CheckNextStep returned a step length that is not kInfinity but the step is indeed boundary limited. Are we actually landing on a boundary?";
	G4Exception("G4CMPBogoliubovQPRandomWalkTransport::AlongStepDoIt", "BogoliubovQPRandomWalkTransport00X",FatalException, msg);
      }
      else{
	G4cout << "---> REL/EY: CheckNextStep shows that we've hit a boundary with our fPathLength computed in G4CMPBogoliubovRandomWalkTransport class. Recalculating time-to-wall and required velocity for this, but not setting new position -- will let Transportation set this." << G4endl;

	/*
	//Compute a representative deltaT at which we reach the wall given diffusion. Note that here, we really should be drawing from a
	//distribution, since there's not a one-to-one for the time at which diffusion from one point to another occurs. But for now we'll
	//make it a single time for simplicity.
	G4double newDiffusionTime = nextStepLength*nextStepLength / 2.0 / fDiffConst;
	*/

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
	fParticleChange.ProposeMomentumDirection(fNewDirection);      
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


/*
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//"Base" version
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
*/
      

  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* 
G4CMPBogoliubovQPRandomWalkTransport::PostStepDoIt(const G4Track& track, const G4Step&)
{
  
  G4cout << "--->REL/EY: In poststepdoit." << G4endl;
  
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
  
  //Pseudocode:
  //1. Determine if we end up close to a boundary using Get2DSafety
  G4cout << "--->REL/EY: In poststepdoit: calculating 2D safety." << G4endl;
  G4double the2DSafety = G4CMP::Get2DSafety(track.GetStep()->GetPreStepPoint()->GetTouchable(),
					    track.GetPosition(),
					    track.GetMomentumDirection(),
					    false);
  G4cout << "--->REL/EY: In poststepdoit: 2D safety calculated to be " << the2DSafety << G4endl;
  
  //2. If we are within epsilon of a boundary, the next step should be made directly into the boundary so that boundary
  //   processes can run. We will therefore not randomize the direction of the next step. Instead, we need to find the direction
  //   toward the boundary in order to figure out how to angle the next vector. Note that an exception to this runs if we're *currently*
  //   on a boundary (i.e. we're leaving from the boundary and the step is at a steep enough angle not to leave the epsilon region around
  //   the surface.
  if( the2DSafety < epsilon ){
    G4cout << "--->REL/EY: In poststepdoit: safety is smaller than epsilon, and finding direction to nearby boundary." << G4endl;
    G4ThreeVector returnDir = FindDirectionToNearbyBoundary(track,the2DSafety);
    fParticleChange.ProposeMomentumDirection(returnDir.unit());
  }
  //3. If we're not within epsilon, then just purely randomize.
  else{
    //Randomize the final state momentum
    G4ThreeVector returnDir = G4RandomDirection();
    returnDir.setZ(0);
    G4cout << "--->REL/EY: In poststepdoit: safety is larger than epsilon, and we're just randomizing the next direction." << G4endl;
    fParticleChange.ProposeMomentumDirection(returnDir.unit());
  }

    
  ClearNumberOfInteractionLengthLeft();		// All processes should do this! 
  return &fParticleChange;

}

//Using the information about the track and the safety computed at the track point, we can find a direction to the
//nearby boundary.
G4ThreeVector G4CMPBogoliubovQPRandomWalkTransport::FindDirectionToNearbyBoundary(const G4Track& track, const G4double the2DSafety ){

  G4cout << "---> REL/EY: Inside FindDirectionToNearbyBoundary" << G4endl;
  
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

  G4cout << "---> REL/EY: pos: " << track.GetPosition() << ", shiftedPoint: " << shiftedPoint << ", original safety: " << the2DSafety << ", shiftedPoint2Dsafety: " << shiftedPoint2DSafety << G4endl;

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
  G4cout << "---> REL/EY: theta for rotation: " << theta << G4endl;

  //Now, we have the theta, but this theta could be in either direction relative to the surface. So what we'll do is try rotating our momDir vector
  //by theta in both directions, stepping the point forward by each, and re-finding safeties. This calculation is specific to the XY plane and
  //should be generalized (REL).
  G4ThreeVector option1 = momDir.rotateZ(theta);
  G4ThreeVector option2 = momDir.rotateZ(-2*theta);
  momDir.rotateZ(theta);
  G4cout << "---> REL/EY: option1: " << option1 << ", option2: " << option2 << G4endl;

  //Recheck safety. Here, we compute the smaller of the two original safeties and make sure we set our
  //"probe" vectors' magnitudes no larger than that. (Here, the 0.9 is to keep it from being exactly one of the safeties)
  G4double smallerSafety = the2DSafety;
  if( shiftedPoint2DSafety < the2DSafety ) smallerSafety = shiftedPoint2DSafety;
  smallerSafety *= 0.9; 
  G4ThreeVector newPosOption1 = track.GetPosition() + smallerSafety*option1;
  G4ThreeVector newPosOption2 = track.GetPosition() + smallerSafety*option2;

  //Now compute the safeties of our two probe vectors' endpoints to see which is closer.
  G4cout << "---> REL/EY: new pos option 1: " << newPosOption1 << G4endl;
  G4cout << "---> REL/EY: new pos option 2: " << newPosOption2 << G4endl; 
  G4double option1Safety = G4CMP::Get2DSafety(track.GetStep()->GetPreStepPoint()->GetTouchable(),
					      newPosOption1,
					      track.GetMomentumDirection(),
					      false);
  G4double option2Safety = G4CMP::Get2DSafety(track.GetStep()->GetPreStepPoint()->GetTouchable(),
					      newPosOption2,
					      track.GetMomentumDirection(),
					      false);

  G4cout << "--->REL/EY: In findDirectionToNearbyBoundary, option 1 safety is " << option1Safety << ", and option 2 safety is: " << option2Safety << G4endl;

  //If option 1 safety is lower, it means we return option 1 as the direction to the boundary
  if( option1Safety < option2Safety ){
    G4cout << "---> REL/EY: Direction to the boundary is option 1: " << option1 << G4endl;
    return option1;
  }
  else if( option2Safety < option1Safety ){
    G4cout << "---> REL/EY: Direction to the boundary is option 2: " << option2 << G4endl;
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
              const G4Track& track, G4double previousStepSize, G4ForceCondition* condition)
{
  G4cout << "REL Here in GetMeanFreePath for the G4CMPBogoliubovQPRandomWalkTransport." << G4endl;

  G4cout << "---> REL/EY MFP: At beginning of MFP, track volume: " << track.GetVolume()->GetName() << G4endl;
  G4cout << "---> REL/EY MFP: At beginning of MFP, previousStepSize: " << previousStepSize << G4endl;
  G4cout << "---> REL/EY MFP: At beginning of MFP, momentum information: " << track.GetMomentumDirection() << G4endl;
  G4cout << "---> REL/EY MFP: At beginning of MFP, position: " << track.GetPosition() << G4endl;
  G4cout << "---> REL/EY MFP: At beginning of MFP, global time: " << track.GetGlobalTime() << G4endl;


  
  //This needs to be done so that we can update the SCUtils information, and since GetMeanFreePath for the discrete bit of this process should (?)
  //run first, it will hopefully happen before the AlongStepGPIL runs and requests the gap.
  if (UpdateMeanFreePathForLatticeChangeover(track)){
    UpdateSCAfterLatticeChange();
  }

  //--------------------------------------------------------------------
  //0. Some preliminaries: setting up safety helper
  if(!fSafetyHelper) {
    G4TransportationManager* transportMgr ;
    transportMgr = G4TransportationManager::GetTransportationManager() ;
    fSafetyHelper = transportMgr->GetSafetyHelper();        
    fSafetyHelper->InitialiseHelper();
  }

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
  G4cout << "---> REL/EY MFP: momentumDir in ASGPIL: " << momentumDir << G4endl;
  G4cout << "---> REL/EY MFP: CurrentVolume by track: " << currentVolume->GetName() << ", CurrentVolumePlusEps, by GetVolumeAtPoint: " << currentVolPlusEps->GetName() << ", currentVolumeMinusEps, by GetVolumeAtPoint: " << currentVolMinEps->GetName() << ",  step status: " << track.GetStep()->GetPreStepPoint()->GetStepStatus() << G4endl;

  //Get energy and velocity of track
  G4double energy = track.GetKineticEnergy();
  G4double velocity = track.GetVelocity();

  //Initialize the lattice manager
  G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();

  
  //If the volume plus epsilon (in the direction of the step) and the current volume are not the same (and if we're
  //in a boundary-limited step), then we're in a turnaround step. This assumes that the real distance to the next
  //geometric feature is not smaller than the boundTolerance -- in that case, this may break down a bit. Return
  //that this is inactive.
  fTrackOnBoundary = false;
  if( currentVolPlusEps != currentVolume && theStatus == fGeomBoundary ){    
    G4cout << "---> REL/EY MFP: In a turnaround step. Killing the transport GPIL." << G4endl;
    fTrackOnBoundary = true;
    isActive = false;
    return DBL_MAX;
  }

  //We can be on a boundary and still be "active" if we're in the volume into which we're moving.
  if( theStatus == fGeomBoundary ){
    G4cout << "---> REL/EY MFP: On a boundary but not in a turnaround step." << G4endl;
    fTrackOnBoundary = true;
  }
  
  G4cout << "---> REL/EY MFP: Not in a turnaround step. Continuing the transport GPIL." << G4endl;

  G4cout << "---> REL/EY MFP: theStatus: " << theStatus << G4endl;
  //--------------------------------------------------------------------  
  //2. Verify that, if we're on a boundary, the direction of momentum is going into a volume with a lattice. Since we've updated the
  //   lattice at this point, we should confirm that the current procUtils lattice is not null, and that it is the
  //   same as the lattice corresponding to what the latticeManager sees as belonging to currentVolPlusEps.
  //   This is purely a check.
  if( theStatus == fGeomBoundary && !LM->HasLattice(currentVolPlusEps) ){ 
    G4ExceptionDescription msg;
    msg << "We're on a boundary (i.e. our boundary is now behind us) and find no lattice in the region where a BogoliubovQP is intending to go during the QP random walk transport's GetMeanFreePath function.";
    G4Exception("G4CMPBogoliubovQPRandomWalkTransport::GetMFP", "BogoliubovQPRandomWalkTransport001",FatalException, msg);
  }

  
  //3. The diffusion constant and other superconductor-specific information are accessed via the
  //   SCUtils class that this inherits from. No need to actually grab it from somewhere: it's already
  //   a data member.
  G4cout << "---> REL/EY MFP: The gap energy, drawn from SCUtils, is: " << fGapEnergy << G4endl;
  G4cout << "---> REL/EY MFP: The Dn, drawn from SCUtils, is: " << fDn << G4endl;
  G4cout << "---> REL/EY MFP: The Teff, drawn from SCUtils, is: " << fTeff << G4endl;  
  G4cout << "---> REL/EY MFP: Gap energy: " << fGapEnergy << ", energy: " << energy << ", DN: " << fDn << G4endl; 
  
  //4. Given this, compute a step length
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
    G4cout << "--->REL/EY: Diffusion constant (adjusted): " << fDiffConst << G4endl;
    
    //Using this diffusion constant, draw a time to reach the boundary
    //established by the safety. For now we use something simple and wrong to test,
    //taking the 2D safety to be the sigma of the diffusion broadening
    //but need the true expression from Grebenkov
    G4double timeStepToBoundary = the2DSafety*the2DSafety / 2.0 / fDiffConst;
    fTimeStepToBoundary = timeStepToBoundary;
    G4cout << "---> REL/EY MFP: Time step to boundary = " << fTimeStepToBoundary << G4endl;
    
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
    G4cout << "---> REL/EY MFP: Setting the number of interaction lengths for RWTransport's Discrete component to 1, and thisMFP = " << thisMFP << G4endl;
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



  /*  What was here first
  *condition = Forced;
  return DBL_MAX;
  */
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

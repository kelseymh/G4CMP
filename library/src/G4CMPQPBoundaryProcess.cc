/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPQPBoundaryProcess.cc
/// \brief Implementation of the G4CMPQPBoundaryProcess class
//
// This process handles the interaction of QPs with
// boundaries.

#include "G4CMPQPBoundaryProcess.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPGeometryUtils.hh"
#include "G4CMPPhononTrackInfo.hh"
#include "G4CMPSurfaceProperty.hh"
#include "G4CMPTrackUtils.hh"
#include "G4CMPVTrackInfo.hh"
#include "G4CMPUtils.hh"
#include "G4ExceptionSeverity.hh"
#include "G4LatticePhysical.hh"
#include "G4ParallelWorldProcess.hh"
#include "G4ParticleChange.hh"
#include "G4PhysicalConstants.hh"
#include "G4LatticeManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include "G4VParticleChange.hh"
#include "G4VSolid.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"



// Constructor
G4CMPQPBoundaryProcess::G4CMPQPBoundaryProcess(const G4String& aName)
  : G4CMPVQPProcess(aName, fQPBoundaryProcess),G4CMPBoundaryUtils(this),procName("G4CMPQPBoundaryProcess") {
}

// Destructor
G4CMPQPBoundaryProcess::~G4CMPQPBoundaryProcess() {
}

// Compute and return step length
G4double G4CMPQPBoundaryProcess::
PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
                                     G4double previousStepSize,
                                     G4ForceCondition* condition) {
  return GetMeanFreePath(aTrack, previousStepSize, condition);
}

G4double G4CMPQPBoundaryProcess::GetMeanFreePath(const G4Track& aTrack,
						 G4double /*prevStepLength*/,
						 G4ForceCondition* condition) {
  //Debugging
  if (verboseLevel > 5) {
    G4cout << "-- G4CMPQPBoundaryProcess::GetMeanFreePath() --" << G4endl;
  }
  
  //Update the lattice so that this process knows about any changes
  UpdateMeanFreePathForLatticeChangeover(aTrack);

  //Use information about the lattice to determine if the QP has
  //been spawned in a self-consistent way. If it's not self-consistent on step
  //1, then throw a fatal exception.
  double eKin = aTrack.GetKineticEnergy();
  G4VPhysicalVolume * currentVol = aTrack.GetVolume();
  G4double stepNumber = aTrack.GetCurrentStepNumber();
  if (stepNumber == 1) {
    if (!IsValidQPVolume(currentVol,eKin)) {
      G4ExceptionDescription msg;
      msg << "Noticed that for the first step, our QP is either not in a "
	  << "superconductor or that the QP energy, " << eKin / eV
	  << " eV, is less than the current volume's gap. You're spawning a "
	  << "quasiparticle either with too low an energy or in the wrong spot "
	  << "for physical accuracy.";
      G4Exception("G4CMPQPBoundaryProcess::GetMeanFreePath",
		  "QPBoundaryProcess001",FatalException, msg);
    }
  }
  
  *condition = Forced;
  return DBL_MAX;
}

//Checks to see if the current volume is a valid one in which a QP with a given
//energy may live. This has three criteria:
//1. Lattice need to exist
//2. Gap0Energy and Tcrit must not be set to their defaults (0)
//3. The kinetic energy of the QP considered (argument 2) must be larger than
//the gap.
G4bool G4CMPQPBoundaryProcess::IsValidQPVolume(G4VPhysicalVolume* volume,
					       G4double qpEKin ) {
  //Debugging
  if (verboseLevel > 5) {
    G4cout << "-- G4CMPQPBoundaryProcess::IsValidQPVolume() --" << G4endl;
  }
  
  
  //Get the lattices from the physical volumes
  //Lattice manager
  G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();
  
  //first lets just check if the volume has a lattice (condition 1)
  if (!LM->HasLattice(volume)) return false;
  
  //Now we need to check to understand if this lattice has existent gap and
  //Tcrit parameters. Philosophically, I think that since we're in the
  //boundary process, I don't want to have/make reference to the parameters
  //needed purely by any of the other processes, like the ElScatMFP, Teff,
  //Tau0_qp, Tau0_ph, or Dn. So here, we check to see that the Tcrit and
  //Gap0Energy are both set to something reasonable. In the absence of other
  //processes which we should be able to turn off, this function should still
  //run. (Condition 2)
  G4LatticePhysical* theLat;
  theLat = LM->GetLattice(volume);
  G4double Gap0Energy = theLat->GetSCDelta0();
  G4double Tcrit = theLat->GetSCTcrit();
  G4double Teff = theLat->GetSCTeff();
  if (Gap0Energy == 0.0 || Tcrit == 0.0) return false;

  //Calculate the nonzero-temperature gap from these using the SCUtils class.
  G4double GapEnergy = ComputeTestGapEnergyAtNonzeroT(Teff,Tcrit,Gap0Energy);

  //Debugging
  if (verboseLevel > 5) {
    G4cout << "RWBoundary IVQPV Function Point A | gapEnergy is "
	   << GapEnergy / CLHEP::eV << " eV" << G4endl;
  }
  
  //Condition 3
  if (GapEnergy > qpEKin) { return false; }
  
  //If we pass all of these, return true
  return true;
}

// Check the pre- and post-step volumes, determine if they have lattices, and
//if they do, set the gaps before and after. If at least one is a valid QP
//volume, then return true. Here, whether the volume is valid also depends on
//the energy of the QP. If the QP energy is smaller than the gap, that region
//should be treated the same as if it's a volume without a crystal. (Are there
//some edge cases that I'm not seeing?)
G4bool G4CMPQPBoundaryProcess::CheckQPVolumes(const G4Step& aStep)
{
  //Check if pre/post step volumes have valid QP volumes
  //Get the lattices from the physical volumes
  //Lattice manager  
  G4double qpEKin = aStep.GetTrack()->GetKineticEnergy();
  preQPVolume =
    IsValidQPVolume(aStep.GetPreStepPoint()->GetPhysicalVolume(),qpEKin);  
  postQPVolume =
    IsValidQPVolume(aStep.GetPostStepPoint()->GetPhysicalVolume(),qpEKin);
  
  //Keep in mind that during turnaround steps, the preQPVolume should be false
  //and the postQPVolume should be true
  return (preQPVolume || postQPVolume);
}

// Process action
G4VParticleChange*
G4CMPQPBoundaryProcess::PostStepDoIt(const G4Track& aTrack,
				     const G4Step& aStep) {

  verboseLevel = G4CMPConfigManager::GetVerboseLevel();

  //Debugging
  if (verboseLevel > 5) {
    G4cout << "-- G4CMPQPBoundaryProcess::PostStepDoIt() --" << G4endl;
    G4cout << "RWBoundary PSDI Function Point A | poststeppoint velocity in "
	   << "RWBoundary poststepdoit is: "
	   << aStep.GetPostStepPoint()->GetVelocity() << G4endl;
    G4cout << "RWBoundary PSDI Function Point A | track velocity in RWBoundary "
	   << " poststepdoit is: " << aTrack.GetVelocity() << G4endl;
  }
  
  // NOTE:  G4VProcess::SetVerboseLevel is not virtual!  Can't overlaod it
  G4CMPBoundaryUtils::SetVerboseLevel(verboseLevel);
  aParticleChange.Initialize(aTrack);

  //Debugging
  if (verboseLevel > 5) {
    G4cout << "RWBoundary PSDI Function Point B | poststeppoint velocity in "
	   << "RWBoundary poststepdoit, after initializing aParticleChange, "
	   << "is: " << aStep.GetPostStepPoint()->GetVelocity() << G4endl;
    G4cout << "RWBoundary PSDI Function Point B | track velocity in RWBoundary "
	   << "poststepdoit, after initializing aParticleChange, is: "
	   << aTrack.GetVelocity() << G4endl;
  }

  //Do a boundary check just as for phonon dynamics
  G4bool checkBoundary = IsGoodBoundary(aStep);
  
  //Debugging
  if (verboseLevel > 5) {
    G4cout << "RWBoundary PSDI Function Point C | After IsGoodBoundary, value "
	   << checkBoundary << G4endl;
  }
  
  //After a boundary check, we also want to do a QP-specific check of the
  //volumes to make sure we understand the relationship between the pre- and
  //post-boundary superconducting gaps. This updates those gap values. First,
  //debugging
  if (verboseLevel > 5) {
    G4cout << "RWBoundary PSDI Function Point D | Before CheckQPVolumes, "
	   << "PostStepDoIt" << G4endl;
  }  
  G4bool checkQPVolumes = CheckQPVolumes(aStep);

  //Debugging
  if (verboseLevel > 5) {
    G4cout << "RWBoundary PSDI Function Point E | After CheckQPVolumes (value "
	   << checkQPVolumes << "), PostStepDoIt" << G4endl;
  }
  if (verboseLevel > 2) {
    G4cout << "G4CMPQPBoundaryProcess: inside PostStepDoIt Check qp volumes "
	   << "result :  " <<checkQPVolumes << G4endl;
    G4cout << "G4CMPQPBoundaryProcess: inside PostStepDoIt Check boundary "
	   << "result :  " <<checkBoundary << G4endl;
  }
  
  //If boundaries or QP volumes aren't satisfied, just return the default
  //post-step do it for discrete processes.
  if (!checkBoundary || !checkQPVolumes) {
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }
  if (verboseLevel>1) {
    G4cout << GetProcessName() << "::PostStepDoIt" << G4endl;
  }

  //Otherwise, apply a boundary action (reflection, absorption, transmission)
  if (verboseLevel > 5) {   
    G4cout << "RWBoundary PSDI Function Poing F | Applying boundary action, "
	   << "PostStepDoIt" << G4endl;
  }  
  ApplyBoundaryAction(aTrack, aStep, aParticleChange);
  ClearNumberOfInteractionLengthLeft();	// All processes should do this!
  
  return &aParticleChange;
}




// Do a reflection depending on the gap conditions between multiple lattices,
//and also depending on the assigned reflection probability given to the
//boundary. Note that here the "gap contition" is now bundled into the
//"postQPVolume" variable. Any volume without a postQPVolume (no lattice,
//default gap, or gap>qpEnergy) will induce reflection.
G4bool G4CMPQPBoundaryProcess::ReflectTrack(const G4Track& /*aTrack*/,
					    const G4Step&) const {
  
  //Note that this is still just blindly copied from the phononboundary action
  // this needs to be customized for QP dynamics ...
  G4double reflProb = GetMaterialProperty("reflProb");
  if (verboseLevel>2) G4cout << " ReflectTrack: reflProb " << reflProb << G4endl;
  //  check if the next volume is a QP lattice if not reflect
  if (!postQPVolume) reflProb =1;
  
  return (G4UniformRand() <= reflProb);
}


// Do absorption of a quasiparticle. The way this is currently run is
// to pass the buck to the BoundaryUtils class, which does a NIEL
// calculation and passage to the partitioner. If we want QPs not to
// go through that process, then we can switch to just the final
// two lines of this function that are commented out.
void G4CMPQPBoundaryProcess::DoAbsorption(const G4Track& aTrack,
					  const G4Step& aStep,
					  G4ParticleChange& aParticleChange) {
  //Debugging
  if (verboseLevel > 5) {
    G4cout << "-- G4CMPQPBoundaryProcess::DoAbsorption() --" << G4endl;
  }
  G4CMPBoundaryUtils::DoAbsorption(aTrack,aStep,aParticleChange);
  
  /* 
  aParticleChange.ProposeTrackStatus(fStopAndKill);
  aParticleChange.ProposeEnergy(0.);
  */
}



// Do reflection of a quasiparticle
void G4CMPQPBoundaryProcess::DoReflection(const G4Track& aTrack,
					  const G4Step& aStep,
					  G4ParticleChange& aParticleChange) {

  //Debugging
  if (verboseLevel > 5) {
    G4cout << "-- G4CMPQPBoundaryProcess::DoReflection() --" << G4endl;
    G4cout << "RWBoundary DR Function Point A | Using reflection where all "
	   << "returned directions are surface norms." << G4endl;
  }
  
  //This function is to be used with QP diffusion. It *will* return the
  //momentum as the surface normal in the direction of motion, as that
  //information is needed by/used by the diffusion class. Check to make sure
  //we're on a volume boundary before attempting reflection.
  G4ThreeVector surfacePoint;
  if (!CheckStepBoundary(aStep, surfacePoint)) {
    if (verboseLevel > 5) {
      G4cout << "RWBoundary DR Function Point B | checking step boundary failed"
	     << " in DoReflection" << G4endl;
    }
    if (verboseLevel>2)
      G4cout << " Boundary point moved to " << surfacePoint << G4endl;
    aParticleChange.ProposePosition(surfacePoint);	// IS THIS CORRECT?!?
  }
  
  if (verboseLevel>1) {
    G4cout << procName << ": Track reflected "
	   << G4CMP::GetTrackInfo<G4CMPVTrackInfo>(aTrack)->ReflectionCount()
	   << " times." << G4endl;
  }
  
  //To determine new random direction, need to know relationship between current
  //direction and the surface normal. If they are more parallel, then we need
  //to ensure that the new direction dotted into the norm is negative. If they
  //are more antiparallel, we need to make sure that the new direction dotted into the norm is positive.    
  G4ThreeVector norm = G4CMP::GetSurfaceNormal(aStep);    
  G4ThreeVector pdir = aTrack.GetMomentumDirection();
  G4ThreeVector newDir;
  
  //your new momentum is very parallel to a surface. Hopefully shouldn't change
  //stuff too much. If initial momentum is in the direction of the surface
  //normal, the return direction should be just the negative of the surface
  //normal
  if (pdir.dot(norm) > 0) {
    newDir = -1*norm;
  } else if (pdir.dot(norm) < 0) {
    newDir = norm;
  } else{
    G4ExceptionDescription msg;
    msg << "Somehow the incoming momentum is exactly parallel to the surface "
	<< "norm? What?";
    G4Exception((GetProcessName()+"::DoReflection").c_str(),
		"G4CMPQPBoundaryProcess002",
		FatalException,
		msg);
  }
  
  //Debugging
  if (verboseLevel > 5) {
    G4cout << "RWBoundary DR Function Point C | inside DoReflection initial "
	   << "direction  " << pdir << G4endl;
    G4cout << "RWBoundary DR Function Point C | inside DoReflection reflected "
	   << "direction  " << newDir << G4endl;
  }
  
  aParticleChange.ProposeMomentumDirection(newDir);
}



// Do transmission of a quasiparticle
void G4CMPQPBoundaryProcess::DoTransmission(const G4Track& aTrack,
					    const G4Step& aStep,
					    G4ParticleChange& aParticleChange) {
  //Debugging
  if (verboseLevel > 5) {
    G4cout << "-- G4CMPQPBoundaryProcess::DoTransmission() --" << G4endl;
  }
  if (verboseLevel > 1) {
    G4cout << procName << ": Track transmission requested" << G4endl;
  }

  //Double-check that you have a proper QP volume in the post-step point. This
  //should never pass, but is a failure mode we should monitor for a bit during
  //debugging.
  if (!postQPVolume) {

    //Debugging
    if (verboseLevel > 5) {
      G4cout << "RWBoundary DT Function Point A | Killing QP inside "
	     << "DoTransmission - postQPVolume is not valid should have been "
	     << "caught in ReflectTrack()!" << G4endl;
    }
    G4ExceptionDescription msg;
    msg << "Noticed that the post-step volume isn't a good QP volume. There is "
	<< "a bug somewhere that needs to be fixed.";
    G4Exception("G4CMPQPBoundaryProcess::DoTransmission",
		"QPBoundaryProcess003",JustWarning, msg);
    DoSimpleKill(aTrack, aStep, aParticleChange);
  }
    
  // Check whether step has proper boundary-stopped geometry
  G4ThreeVector surfacePoint;
  if (!CheckStepBoundary(aStep, surfacePoint)) {

    //Debugging
    if (verboseLevel > 5) {
      G4cout << "RWBoundary DT Function Point B | REL checking step boundary "
	     << "failed in DoTransmission" << G4endl;
    }
    if (verboseLevel>2)
      G4cout << " Boundary point moved to " << surfacePoint << G4endl;
    
    aParticleChange.ProposePosition(surfacePoint);    // IS THIS CORRECT?!?
  }
  
  //Since the lattice hasn't changed yet, change it here. (This also happens
  //at the MFP calc point at the beginning of the next step, but it's nice to
  //have it here so we can use the new lattice info to help figure out vdir,
  //etc.)
  this->SetLattice(G4LatticeManager::GetLatticeManager()->GetLattice(aStep.GetPostStepPoint()->GetPhysicalVolume()));
  UpdateSCAfterLatticeChange();
    
  G4ThreeVector vdir = aTrack.GetMomentumDirection();
  G4ThreeVector norm = G4CMP::GetSurfaceNormal(aStep);

  //Make sure that the norm and the vdir are in the same direction so we return
  //a norm that is in the direction of travel
  if (vdir.dot(norm) < 0) norm = -1*norm;
  
  G4int transmissionType = 1;
  
  //For nice visualization/debugging
  if (transmissionType == 0) { //REL THIS SHOULD BE DEPRECATED  
    aParticleChange.ProposeMomentumDirection(vdir);
  } else if (transmissionType == 1) { //REL THIS SHOULD BE THE ONLY ONE HERE.
    //^For use with QP transport
    aParticleChange.ProposeMomentumDirection(norm);
  } else { //Huh?
    G4ExceptionDescription msg;
    msg << "Noticed that we are using an undefined transmissionType in "
	<< "DoTransmission for QPs. Fix.";
    G4Exception("G4CMPQPBoundaryProcess::DoTransmission",
		"QPBoundaryProcess004",FatalException, msg);
  }
}




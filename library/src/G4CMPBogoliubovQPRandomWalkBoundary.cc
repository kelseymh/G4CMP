/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPBogoliubovQPRandomWalkBoundary.cc
/// \brief Implementation of the G4CMPBogoliubovQPRandomWalkBoundary class
//
// This process handles the interaction of BogoliubovQPs with
// boundaries.

#include "G4CMPBogoliubovQPRandomWalkBoundary.hh"
#include "G4CMPSCUtils.hh"
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
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VParticleChange.hh"
#include "G4VSolid.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Constructor
G4CMPBogoliubovQPRandomWalkBoundary::G4CMPBogoliubovQPRandomWalkBoundary(const G4String& aName)
  : G4VBogoliubovQPProcess(aName, fBogoliubovQPRandomWalkBoundary),G4CMPBoundaryUtils(this),procName("G4CMPBogloliubovQPRandomWalkBoundary")
{
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Destructor
G4CMPBogoliubovQPRandomWalkBoundary::~G4CMPBogoliubovQPRandomWalkBoundary()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Compute and return step length
G4double G4CMPBogoliubovQPRandomWalkBoundary::
PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
                                     G4double previousStepSize,
                                     G4ForceCondition* condition)
{
  return GetMeanFreePath(aTrack, previousStepSize, condition);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4double G4CMPBogoliubovQPRandomWalkBoundary::GetMeanFreePath(const G4Track& aTrack,
							      G4double /*prevStepLength*/,
							      G4ForceCondition* condition)
{
  G4cout << "REL -- G4CMPQPRandomWalkBoundary::GetMeanFreePath()" << G4endl;
  
  //Update the lattice so that this process knows about any changes
  UpdateMeanFreePathForLatticeChangeover(aTrack);

  //Use information about the lattice to determine if the BogoliubovQP has been spawned in a self-consistent way. If it's not
  //self-consistent on step 1, then throw a fatal exception.
  double eKin = aTrack.GetKineticEnergy();
  G4VPhysicalVolume * currentVol = aTrack.GetVolume();
  G4double stepNumber = aTrack.GetCurrentStepNumber();
  if( stepNumber == 1 ){
    if( !IsValidQPVolume(currentVol,eKin) ){
      G4ExceptionDescription msg;
      msg << "Noticed that for the first step, our QP is either not in a superconductor or that the QP energy, " << eKin / eV << " eV, is less than the current volume's gap. You're spawning a quasiparticle either with too low an energy or in the wrong spot for physical accuracy.";
      G4Exception("G4CMPBogoliubovQPRandomWalkBoundary::GetMeanFreePath", "BogoliubovQPRandomWalkBoundary001",FatalException, msg);
    }
  }
  
  *condition = Forced;
  return DBL_MAX;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//Checks to see if the current volume is a valid one in which a QP with a given energy may live. This has three criteria:
//1. Lattice need to exist
//2. Gap0Energy and Tcrit must not be set to their defaults (0)
//3. The kinetic energy of the QP considered (argument 2) must be larger than the gap.
//Question for Eric: is there ever a time when we'll want to pull condition 3 OUT of the statement about whether a QP
//can exist in a volume? 
G4bool G4CMPBogoliubovQPRandomWalkBoundary::IsValidQPVolume(G4VPhysicalVolume* volume, G4double qpEKin )
{
  //Get the lattices from the physical volumes
  //Lattice manager
  G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();
  
  //first lets just check if the volume has a lattice (condition 1)
  if (!LM->HasLattice(volume)){ return false; }
  
  //Now we need to check to understand if this lattice has existent gap and Tcrit parameters.
  //Philosophically, I think that since we're in the boundary process, I don't want to have/make
  //reference to the parameters needed purely by any of the other processes, like the ElScatMFP,
  //Teff, Tau0_qp, Tau0_ph, or Dn. So here, we check to see that the Tcrit and Gap0Energy are both set
  //to something reasonable. In the absence of other processes which we should be able to turn off,
  //this function should still run. (Condition 2)
  G4LatticePhysical* theLat;
  theLat = LM->GetLattice(volume);
  G4double Gap0Energy = theLat->GetSCDelta0();
  G4double Tcrit = theLat->GetSCTcrit();
  G4double Teff = theLat->GetSCTeff();
  if( Gap0Energy == 0.0 || Tcrit == 0.0 ) return false;

  //Calculate the nonzero-temperature gap from these using the SCUtils class.
  G4double GapEnergy = ComputeTestGapEnergyAtNonzeroT(Teff,Tcrit,Gap0Energy);
  G4cout << "REL -- in G4CMPBogoliubovQPRandomWalkBoundary::IsValidQPVolume(): gapEnergy is " << GapEnergy / CLHEP::eV << " eV" << G4endl;
  
  //Condition 3
  if( GapEnergy > qpEKin ){ return false; }
  
  //If we pass all of these, return true
  return true;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Check the pre- and post-step volumes, determine if they have lattices, and if they
// do, set the gaps before and after. If at least one is a valid QP volume, then
// return true. Here, whether the volume is valid also depends on the energy of the QP.
// If the QP energy is smaller than the gap, that region should be treated the same as if
// it's a volume without a crystal. (Are there some edge cases that I'm not seeing?)
G4bool G4CMPBogoliubovQPRandomWalkBoundary::CheckQPVolumes(const G4Step& aStep)
{
  //Check if pre/post step volumes have valid QP volumes
  //Get the lattices from the physical volumes
  //Lattice manager
  G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();
  G4double qpEKin = aStep.GetTrack()->GetKineticEnergy();
  preQPVolume = IsValidQPVolume(aStep.GetPreStepPoint()->GetPhysicalVolume(),qpEKin);  
  postQPVolume = IsValidQPVolume(aStep.GetPostStepPoint()->GetPhysicalVolume(),qpEKin);
  
  //Keep in mind that during turnaround steps, the preQPVolume should be false and the postQPVolume should be true
  return (preQPVolume || postQPVolume);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Process action
G4VParticleChange*
G4CMPBogoliubovQPRandomWalkBoundary::PostStepDoIt(const G4Track& aTrack,
						  const G4Step& aStep) {

  G4cout << "REL -- G4CMPQPRandomWalkBoundary::PostStepDoIt()" << G4endl;
  
  //Note that this is still just blindly copied from the phononboundary action
  // this needs to be customized for QP dynamics ...
  // NOTE:  G4VProcess::SetVerboseLevel is not virtual!  Can't overlaod it
  G4CMPBoundaryUtils::SetVerboseLevel(verboseLevel);
  aParticleChange.Initialize(aTrack);

  //Do a boundary check just as for phonon dynamics
  G4bool checkBoundary = IsGoodBoundary(aStep);
  G4cout << "REL -- After IsGoodBoundary, value " << checkBoundary << G4endl;

  //After a boundary check, we also want to do a QP-specific check of the volumes to make sure we understand the
  //relationship between the pre- and post-boundary superconducting gaps. This updates those gap values.
  G4cout << "REL -- Before CheckQPVolumes, PostStepDoIt" << G4endl;
  G4bool checkQPVolumes = CheckQPVolumes(aStep);
  G4cout << "REL -- After CheckQPVolumes (value " << checkQPVolumes << "), PostStepDoIt" << G4endl;
  if (verboseLevel>2){
    G4cout << "G4CMPBogoliubovQPRandomWalkBoundary: inside PostStepDoIt Check qp volumes result :  " <<checkQPVolumes << G4endl;
    G4cout << "G4CMPBogoliubovQPRandomWalkBoundary: inside PostStepDoIt Check boundary result :  " <<checkBoundary << G4endl;
  }
  
  //If boundaries or QP volumes aren't satisfied, just return the default post-step do it for discrete processes.
  if (!checkBoundary || !checkQPVolumes) return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);  
  if (verboseLevel>1) G4cout << GetProcessName() << "::PostStepDoIt" << G4endl;

  //Otherwise, apply a boundary action (reflection, absorption, transmission)
  G4cout << "REL -- Applying boundary action, PostStepDoIt" << G4endl;
  ApplyBoundaryAction(aTrack, aStep, aParticleChange);
  ClearNumberOfInteractionLengthLeft();		// All processes should do this!
  return &aParticleChange;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Do a reflection depending on the gap conditions between multiple lattices, and
// also depending on the assigned reflection probability given to the boundary. Note that
// here the "gap contition" is now bundled into the "postQPVolume" variable. Any volume without
// a postQPVolume (no lattice, default gap, or gap>qpEnergy) will induce reflection.
G4bool G4CMPBogoliubovQPRandomWalkBoundary::ReflectTrack(const G4Track& aTrack, const G4Step&) const {
  //Note that this is still just blindly copied from the phononboundary action
  // this needs to be customized for QP dynamics ...
  G4double reflProb = GetMaterialProperty("reflProb");
  if (verboseLevel>2) G4cout << " ReflectTrack: reflProb " << reflProb << G4endl;
  //  check if the next volume is a QP lattice if not reflect
  if (!postQPVolume) reflProb =1;
  //    G4cout << "G4CMPBogoliubovQPRandomWalkBoundary: inside ReflectTack! reflection prob :  " <<reflProb << G4endl;
  // Check superconducting gap of the next volume compared to QP energy
  //G4cout << "REL ReflectTrack: postQPVolume " << postQPVolume << G4endl;
  //G4cout << "REL  ReflectTrack: reflProb " << reflProb << G4endl;
  //G4cout << "REL ReflectTrack: postSCGap: " << postSCGap << ", qp energy: " << Eqp << G4endl;
  return (G4UniformRand() <= reflProb);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Do absorption of a quasiparticle.
void G4CMPBogoliubovQPRandomWalkBoundary::DoAbsorption(const G4Track& aTrack,
						       const G4Step&,
						       G4ParticleChange& aParticleChange) {
  //Note that this is still just blindly copied from the G4CMPBoundaryUtils
  // this needs to be customized for QP dynamics ...
  if (verboseLevel>1) G4cout << procName << ": Track absorbed" << G4endl;
  G4double ekin = procUtils->GetKineticEnergy(aTrack);
  G4cout << "REL: in doabsorption, ekin: " << ekin << " for QP." << G4endl;
  aParticleChange.ProposeNonIonizingEnergyDeposit(ekin);
  aParticleChange.ProposeTrackStatus(fStopAndKill);
  aParticleChange.ProposeEnergy(0.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Generate Lambertian/diffuse reflection off internal surfaces
G4ThreeVector G4CMPBogoliubovQPRandomWalkBoundary::GetLambertianVector(const G4ThreeVector& surfNorm) const {
  G4ThreeVector vdir;
  const G4int maxTries = 1000;
  G4int nTries = 0;
  do {
    vdir = G4CMP::LambertReflection(surfNorm);
  } while (nTries++ < maxTries &&
	   (vdir.dot(surfNorm) > 0.0) );
  return vdir;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Do diffuse reflection of a quasiparticle
void G4CMPBogoliubovQPRandomWalkBoundary::DoReflection(const G4Track& aTrack,
						       const G4Step& aStep,
						       G4ParticleChange& aParticleChange)
{
  G4cout << "REL -- G4CMPQPRandomWalkBoundary::DoReflection()" << G4endl;
  
  //REL currently hardcoded but should fix
  bool isLambertian = false;

  if( isLambertian ){

    G4cout << "REL -- G4CMPQPRandomWalkBoundary::Lambertian Reflection()" << G4endl;
    
    // Check whether step has proper boundary-stopped geometry
    G4ThreeVector surfacePoint;
    if (!CheckStepBoundary(aStep, surfacePoint)) {
      G4cout << "REL checking step boundary failed in DoReflection" << G4endl;
      if (verboseLevel>2)
	G4cout << " Boundary point moved to " << surfacePoint << G4endl;
      aParticleChange.ProposePosition(surfacePoint);	// IS THIS CORRECT?!?
    }

    //Get a lambertian/diffuse reflection
    G4ThreeVector surfNorm = G4CMP::GetSurfaceNormal(aStep);
    G4ThreeVector vdir = GetLambertianVector(surfNorm);

    // If reflection failed, report problem and kill the track
    if (vdir.dot(surfNorm) > 0.0 ){
      G4Exception((GetProcessName()+"::DoReflection").c_str(), "G4CMPBogoliubovQPRandomWalkBoundary003",
		  JustWarning, "BogoliubovQP reflection failed");
      DoSimpleKill(aTrack, aStep, aParticleChange);
      return;
    }

    // SANITY CHECK:  Project a 1 um step in the new direction, see if it
    // is still in the correct (pre-step) volume.

    if (verboseLevel>2) {
      G4ThreeVector stepPos = surfacePoint + 1*um * vdir;

      G4cout << " New travel direction " << vdir
	     << "\n from " << surfacePoint << "\n   to " << stepPos << G4endl;

      G4ThreeVector stepLocal = GetLocalPosition(stepPos);
      G4VSolid* solid = aStep.GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetSolid();

      EInside place = solid->Inside(stepLocal);
      G4cout << " After trial step, " << (place==kInside ? "inside"
					  : place==kOutside ? "OUTSIDE"
					  : "on surface") << G4endl;
    }

    aParticleChange.ProposeMomentumDirection(vdir);
  }

  //Specular
  else{

    G4cout << "REL -- G4CMPQPRandomWalkBoundary::Specular Reflection()" << G4endl;

    //Check to make sure we're on a volume boundary before attempting reflection.
    G4ThreeVector surfacePoint;
    if (!CheckStepBoundary(aStep, surfacePoint)) {
      G4cout << "REL checking step boundary failed in DoReflection" << G4endl;
      if (verboseLevel>2)
	G4cout << " Boundary point moved to " << surfacePoint << G4endl;
      aParticleChange.ProposePosition(surfacePoint);	// IS THIS CORRECT?!?
    }
    
    if (verboseLevel>1) {
      G4cout << procName << ": Track reflected "
	     << G4CMP::GetTrackInfo<G4CMPVTrackInfo>(aTrack)->ReflectionCount()
	     << " times." << G4endl;
    }
  
    G4ThreeVector pdir = aTrack.GetMomentumDirection();
    //    G4ThreeVector pdir = aStep.GetPreStepPoint()->GetMomentumDirection();
    //    if (verboseLevel>2)
    G4cout << "G4CMPBogoliubovQPRandomWalkBoundary: inside DoReflection initial direction  " <<pdir << G4endl;    
    
    G4ThreeVector norm = G4CMP::GetSurfaceNormal(aStep);    // Outward normal
    G4cout << "REL norm is " << norm.x() << ", " << norm.y() << ", " << norm.z() << G4endl;
    pdir -= 2.*(pdir.dot(norm))*norm;            // Reverse along normal

    //if (verboseLevel>2)
    G4cout << "G4CMPBogoliubovQPRandomWalkBoundary: inside DoReflection reflected direction  " <<pdir << G4endl;
    
    aParticleChange.ProposeMomentumDirection(pdir);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Do transmission of a quasiparticle
void G4CMPBogoliubovQPRandomWalkBoundary::DoTransmission(const G4Track& aTrack,
							 const G4Step& aStep,
							 G4ParticleChange& aParticleChange) {

  if (verboseLevel>1){
    G4cout << procName << ": Track transmission requested" << G4endl;
  }

  //Double-check that you have a proper QP volume in the post-step point. This should never pass, but is a failure mode we should monitor for
  //a bit during debugging.
  if (!postQPVolume){ 
    G4cout << "Killing QP inside DoTransmission - postQPVolume is not valid should have been caught in ReflectTrack()!" << G4endl;
    G4ExceptionDescription msg;
    msg << "Noticed that the post-step volume isn't a good QP volume. There is a bug somewhere that needs to be fixed.";
    G4Exception("G4CMPBogoliubovQPRandomWalkBoundary::DoTransmission", "BogoliubovQPRandomWalkBoundary004",JustWarning, msg);
    DoSimpleKill(aTrack, aStep, aParticleChange);
  }
    
  // Check whether step has proper boundary-stopped geometry
  G4ThreeVector surfacePoint;
  if (!CheckStepBoundary(aStep, surfacePoint)) {
    G4cout << "REL checking step boundary failed in DoTransmission" << G4endl;
    if (verboseLevel>2)
      G4cout << " Boundary point moved to " << surfacePoint << G4endl;
    
    aParticleChange.ProposePosition(surfacePoint);    // IS THIS CORRECT?!?
  }

  //Since the lattice hasn't changed yet, change it here. (This also happens at the MFP calc point at the beginning of the next step,
  //but it's nice to have it here so we can use the new lattice info to help figure out vdir, etc.)
  this->SetLattice(G4LatticeManager::GetLatticeManager()->GetLattice(aStep.GetPostStepPoint()->GetPhysicalVolume()));
  UpdateSCAfterLatticeChange();
    
  G4ThreeVector vdir = aTrack.GetMomentumDirection();
  aParticleChange.ProposeMomentumDirection(vdir);
}

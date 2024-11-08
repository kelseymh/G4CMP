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

// Constructor and destructor
G4CMPBogoliubovQPRandomWalkBoundary::G4CMPBogoliubovQPRandomWalkBoundary(const G4String& aName)
  : G4CMPVProcess(aName, fBogoliubovQPRandomWalkBoundary),G4CMPBoundaryUtils(this),procName("G4CMPBogloliubovQPRandomWalkBoundary") {;}

G4CMPBogoliubovQPRandomWalkBoundary::~G4CMPBogoliubovQPRandomWalkBoundary() {;}

// Compute and return step length
G4double G4CMPBogoliubovQPRandomWalkBoundary::
PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
                                     G4double previousStepSize,
                                     G4ForceCondition* condition) {
  return GetMeanFreePath(aTrack, previousStepSize, condition);
}

G4double G4CMPBogoliubovQPRandomWalkBoundary::GetMeanFreePath(const G4Track& /*aTrack*/,
                                             G4double /*prevStepLength*/,
                                             G4ForceCondition* condition) {
  *condition = Forced;
  return DBL_MAX;
}

G4bool G4CMPBogoliubovQPRandomWalkBoundary::IsValidQPVolume(G4VPhysicalVolume* volume) {
    //Get the lattices from the physical volumes
    //Lattice manager
    G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();
    //first lets just check if the volume has a lattice
    if (LM->HasLattice(volume)){
        //Now check if it the lattice has a nonzero Gap and Diffusion constant
        G4LatticePhysical* theLat;
        theLat = LM->GetLattice(volume);
        G4double Gap0Energy = theLat->GetSCDelta0();
        G4double Dn = theLat->GetSCDn();
        if ((Dn!=0) && (Gap0Energy!=0)){
            //Then this is a valid QP lattice
            return true;
        }else{
            return false;
        }
    }else{
        //Volume does not have a lattice -> is not a valid QP volume
        return false;
    }
}
    
G4bool G4CMPBogoliubovQPRandomWalkBoundary::CheckQPVolumes(const G4Step& aStep){
    preSCGap = DBL_MAX;
    postSCGap = DBL_MAX;
    //Check if pre/post step volumes have valid QP volumes
    //Get the lattices from the physical volumes
    //Lattice manager
    G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();
    preQPVolume = IsValidQPVolume(aStep.GetPreStepPoint()->GetPhysicalVolume());
    if (preQPVolume){
        //Gap is not stored in in a matieral properties table but in the lattice itself
        G4LatticePhysical* preLattice = LM->GetLattice(aStep.GetPreStepPoint()->GetPhysicalVolume());
        preSCGap = preLattice->GetSCDelta0();
    }

    postQPVolume = IsValidQPVolume(aStep.GetPostStepPoint()->GetPhysicalVolume());
    if (postQPVolume){
        //Gap is not stored in in a matieral properties table but in the lattice itself
        G4LatticePhysical* postLattice = LM->GetLattice(aStep.GetPostStepPoint()->GetPhysicalVolume());
        postSCGap = postLattice->GetSCDelta0();
      }
    return (preQPVolume || postQPVolume);
}

// Process action
G4VParticleChange*
G4CMPBogoliubovQPRandomWalkBoundary::PostStepDoIt(const G4Track& aTrack,
                                         const G4Step& aStep) {
  //Note that this is still just blindly copied from the phononboundary action
  // this needs to be customized for QP dynamics ...
  // NOTE:  G4VProcess::SetVerboseLevel is not virtual!  Can't overlaod it
  G4CMPBoundaryUtils::SetVerboseLevel(verboseLevel);

  aParticleChange.Initialize(aTrack);
  
  G4bool checkBoundary = IsGoodBoundary(aStep);
    
  //the following also updates preSCGap and postSCGap parameters of the class
  G4bool checkQPVolumes = CheckQPVolumes(aStep);
  if (verboseLevel>2){
        G4cout << "G4CMPBogoliubovQPRandomWalkBoundary: inside PostStepDoIt Check qp volumes result :  " <<checkQPVolumes << G4endl;
        G4cout << "G4CMPBogoliubovQPRandomWalkBoundary: inside PostStepDoIt Check boundary result :  " <<checkBoundary << G4endl;
  }
    
  if (!checkBoundary || !checkQPVolumes) return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);

  if (verboseLevel>1) G4cout << GetProcessName() << "::PostStepDoIt" << G4endl;

  ApplyBoundaryAction(aTrack, aStep, aParticleChange);

  ClearNumberOfInteractionLengthLeft();		// All processes should do this!
  return &aParticleChange;
}

G4bool G4CMPBogoliubovQPRandomWalkBoundary::ReflectTrack(const G4Track& aTrack, const G4Step&) const {
    //Note that this is still just blindly copied from the phononboundary action
    // this needs to be customized for QP dynamics ...
  G4double reflProb = GetMaterialProperty("reflProb");
  if (verboseLevel>2) G4cout << " ReflectTrack: reflProb " << reflProb << G4endl;
//  check if the next volume is a QP lattice if not reflect
   if (!postQPVolume) reflProb =1;
//    G4cout << "G4CMPBogoliubovQPRandomWalkBoundary: inside ReflectTack! reflection prob :  " <<reflProb << G4endl;
   // Check superconducting gap of the next volume compared to QP energy
   G4double Eqp = procUtils->GetKineticEnergy(aTrack);
   if (Eqp<postSCGap) reflProb = 1;
  return (G4UniformRand() <= reflProb);
}

void G4CMPBogoliubovQPRandomWalkBoundary::DoAbsorption(const G4Track& aTrack,
                      const G4Step&,
                      G4ParticleChange& aParticleChange) {
    //Note that this is still just blindly copied from the G4CMPBoundaryUtils
    // this needs to be customized for QP dynamics ...
  if (verboseLevel>1) G4cout << procName << ": Track absorbed" << G4endl;
  G4double ekin = procUtils->GetKineticEnergy(aTrack);
  aParticleChange.ProposeNonIonizingEnergyDeposit(ekin);
  aParticleChange.ProposeTrackStatus(fStopAndKill);
  aParticleChange.ProposeEnergy(0.);
}

void G4CMPBogoliubovQPRandomWalkBoundary::DoReflection(const G4Track& aTrack,
                      const G4Step& aStep,
                      G4ParticleChange& aParticleChange) {
  if (verboseLevel>1) {
    G4cout << procName << ": Track reflected "
           << G4CMP::GetTrackInfo<G4CMPVTrackInfo>(aTrack)->ReflectionCount()
       << " times." << G4endl;
  }
  
  G4ThreeVector pdir = aTrack.GetMomentumDirection();
//    G4ThreeVector pdir = aStep.GetPreStepPoint()->GetMomentumDirection();
  if (verboseLevel>2) G4cout << "G4CMPBogoliubovQPRandomWalkBoundary: inside DoReflection initial direction  " <<pdir << G4endl;
  G4ThreeVector norm = G4CMP::GetSurfaceNormal(aStep);    // Outward normal
  pdir -= 2.*(pdir.dot(norm))*norm;            // Reverse along normal

  if (verboseLevel>2)G4cout << "G4CMPBogoliubovQPRandomWalkBoundary: inside DoReflection reflected direction  " <<pdir << G4endl;
    
  aParticleChange.ProposeMomentumDirection(pdir);
}

void G4CMPBogoliubovQPRandomWalkBoundary::DoTransmission(const G4Track& aTrack,
                   const G4Step& aStep,
                   G4ParticleChange& aParticleChange) {

  if (verboseLevel>1){
    G4cout << procName << ": Track transmission requested" << G4endl;
  }
  if (postQPVolume){
      G4cout << "Killing QP inside DoTransmission - postQPVolume is not valid should have been caught in ReflectTrack()!" << G4endl;
      DoSimpleKill(aTrack, aStep, aParticleChange);
  }
  // Check superconducting gap of the next volume compared to QP energy
  G4double Eqp = procUtils->GetKineticEnergy(aTrack);
  if (Eqp>=postSCGap){
      // Check whether step has proper boundary-stopped geometry
      G4ThreeVector surfacePoint;
      if (!CheckStepBoundary(aStep, surfacePoint)) {
          G4cout << "REL checking step boundary failed in DoTransmission" << G4endl;
          if (verboseLevel>2)
            G4cout << " Boundary point moved to " << surfacePoint << G4endl;

          aParticleChange.ProposePosition(surfacePoint);    // IS THIS CORRECT?!?
      }
      G4ThreeVector vdir = aTrack.GetMomentumDirection();
      aParticleChange.ProposeMomentumDirection(vdir);
  }else{
      G4cout << "Killing QP inside DoTransmission - QP energy not enough energy to transport into volume should have been caught in ReflectTrack()" << G4endl;
      DoSimpleKill(aTrack, aStep, aParticleChange);
  }
  
}

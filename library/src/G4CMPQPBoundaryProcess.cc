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
#include "G4CMPQPDiffusionTimeStepperProcess.hh"
#include "G4CMPQPDiffusionTimeStepperRate.hh"
#include "G4CMPQPLocalTrappingProcess.hh"
#include "G4CMPQPLocalTrappingRate.hh"
#include "G4CMPQPRadiatesPhononProcess.hh"
#include "G4CMPQPRadiatesPhononRate.hh"
#include "G4CMPQPRecombinationProcess.hh"
#include "G4CMPQPRecombinationRate.hh"
#include "G4CMPSurfaceProperty.hh"
#include "G4CMPVScatteringRate.hh"
#include "G4CMPTrackUtils.hh"
#include "G4CMPVTrackInfo.hh"
#include "G4CMPUtils.hh"
#include "G4ExceptionSeverity.hh"
#include "G4LatticePhysical.hh"
#include "G4ParallelWorldProcess.hh"
#include "G4ParticleChange.hh"
#include "G4PhysicalConstants.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
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




// Do a reflection based on the gap conditions between multiple lattices, and
// based on a reflection probability rooted in an approximate technique that
// we're applying to un-bias the boundary behavior induced by WoS.
// Note that here the "gap contition" is now bundled into the
//"postQPVolume" variable. Any volume without a postQPVolume (no lattice,
// default gap, or gap>qpEnergy) will induce reflection.
// Note: this function is now updated to CHANGE the lattice to the "far side"
// volume lattice -- this is okay since we're either going to transmit fully
// or transmit into a zero-length step, which will update the lattice to that
// anyway (temporarily)
G4bool G4CMPQPBoundaryProcess::ReflectTrack(const G4Track& aTrack,
                                            const G4Step& aStep) {

  //Take human-installed reflProb
  G4double reflProb = GetMaterialProperty("reflProb");
  if (verboseLevel>2) G4cout << " ReflectTrack: reflProb " << reflProb
                             << G4endl;
  
  //  check if the next volume is a QP lattice if not reflect
  if (!postQPVolume) { return (G4UniformRand() <= 1); }

  //OVERARCHING NARRATIVE OF THIS GEOMETRY-AND-INTERFACE-ADAPTIVE REFLECTION
  //Now, we are going to attempt to undo the bias at transparent interfaces
  //caused by the walk-on-spheres algorithm. The prescription is as follows.
  //1. Assuming we have valid QP volumes in both the pre- and post-step
  //   volumes, we can compute geometrical safeties in both of these. Do so.
  //2. Also compute "effective diffusion distances" for the other potentially-
  //   step-limiting processes on BOTH sides of the interface.
  //3. On each side of the interface, take the shortest of the geometrical
  //   safety and the other diffusion-limiting processes, and compute a
  //   "generalized safety" from the winner: if the geo safety, take it
  //   straight, if something else, use diffusion constant to compute an avg
  //   displacement.
  //   a.) Note that this isn't perfect for the non-geo processes. The "right"
  //       way to do this is to get access to the "number of taus" left for this
  //       process in the process race, but it's hard/janky enough just to get
  //       the (energy dependent) average tau for a process while in a
  //       completely different process class, so I'm going to accept however
  //       this biases us.
  //4. For each generalized safety, estimate the number of *microphysical*
  //   returns to the boundary (on both sides). This is not particularly
  //   intuitive but there is documentation that motivates an approximate
  //   analytical form for this on the G4CMP Confluence (as of a future date).
  //5. Using the relative number of returns-to-interface on both sides of the
  //   interface, we compute an adaptive p_abs, which we weight with the
  //   user's parameter to compute the reflection probability.

  //Start by generating useful information from both sides: position, surface
  //normals, and lattices  
  G4ThreeVector pos = aStep.GetPostStepPoint()->GetPosition();
  G4ThreeVector post_vol_norm =
    aStep.GetPostStepPoint()->GetMomentumDirection(); //Outward
  G4ThreeVector pre_vol_norm = post_vol_norm * -1; //Inward
  const G4LatticePhysical * postLattice = G4LatticeManager::GetLatticeManager()
    ->GetLattice(aStep.GetPostStepPoint()->GetPhysicalVolume());
  const G4LatticePhysical * preLattice = G4LatticeManager::GetLatticeManager()
    ->GetLattice(aStep.GetPreStepPoint()->GetPhysicalVolume());
  
  //Generate the geometrical safeties on both sides of the boundary
  double preStepSafety =
    G4CMP::Get2DSafety(aStep.GetPreStepPoint()->GetTouchable(), pos,
                       aStep.GetPreStepPoint()->GetMomentumDirection(),
                       true, true, pre_vol_norm);
  double postStepSafety =
    G4CMP::Get2DSafety(aStep.GetPostStepPoint()->GetTouchable(), pos,
                       aStep.GetPreStepPoint()->GetMomentumDirection(),
                       true, true, post_vol_norm);

  //Get the other active processes' mean free paths.

  //First, grab the process vector for the current QP and pick out the relevant
  //processes we care about.
  G4ProcessVector pl =
    *(aTrack.GetParticleDefinition()->GetProcessManager()->GetProcessList());

  //Second, check to see which processes are active. Here the convention is
  //0 = qpRadiatesPhonon
  //1 = qpRecombination
  //2 = qpDiffusionTimeStepper
  //3 = qpLocalTrapping
  const int nRelevProcs = 4;
  bool relevantQPProcessesActive[nRelevProcs] = {false,false,false,false};
  int relevantQPProcessIDs[nRelevProcs] = {-1,-1,-1,-1};
  for (int iL = 0; iL < pl.size(); ++iL) {
    if (pl[iL]->GetProcessName() == "qpRadiatesPhonon") {
      relevantQPProcessesActive[0] = aTrack.GetParticleDefinition()
        ->GetProcessManager()->GetProcessActivation(pl[iL]);
      relevantQPProcessIDs[0] = iL;
    }    
    if (pl[iL]->GetProcessName() == "qpRecombination") {
      relevantQPProcessesActive[1] = aTrack.GetParticleDefinition()
        ->GetProcessManager()->GetProcessActivation(pl[iL]);
      relevantQPProcessIDs[1] = iL;
    }
    if (pl[iL]->GetProcessName() == "qpDiffusionTimeStepper") {
      relevantQPProcessesActive[2] = aTrack.GetParticleDefinition()
        ->GetProcessManager()->GetProcessActivation(pl[iL]);
      relevantQPProcessIDs[2] = iL;
    }
    if (pl[iL]->GetProcessName() == "qpLocalTrapping") {
      relevantQPProcessesActive[3] = aTrack.GetParticleDefinition()
        ->GetProcessManager()->GetProcessActivation(pl[iL]);
      relevantQPProcessIDs[3] = iL;
    }
  }

  //Now for active processes, loop through and compute two quantities: the
  //pre-step rates, followed by an update to the lattice, then the post-step
  //rates.
  G4double preVolRates[nRelevProcs] = {-1,-1,-1,-1};
  if (relevantQPProcessesActive[0]) {
    preVolRates[0] = dynamic_cast<G4CMPQPRadiatesPhononProcess*>
      (pl[relevantQPProcessIDs[0]])->GetRateModel()->Rate(aTrack);
  }
  if (relevantQPProcessesActive[1]) {
    preVolRates[1] = dynamic_cast<G4CMPQPRecombinationProcess*>
      (pl[relevantQPProcessIDs[1]])->GetRateModel()->Rate(aTrack);
  }
  if (relevantQPProcessesActive[2]) {
    preVolRates[2] = dynamic_cast<G4CMPQPDiffusionTimeStepperProcess*>
      (pl[relevantQPProcessIDs[2]])->GetRateModel()->Rate(aTrack);
  }
  if (relevantQPProcessesActive[3]) {
    preVolRates[3] = dynamic_cast<G4CMPQPLocalTrappingProcess*>
      (pl[relevantQPProcessIDs[3]])->GetRateModel()->Rate(aTrack);
  }

  
  //Do the transition to the new lattice (used to be in DoTransmission)
  //Since the lattice hasn't changed yet, change it here. (This also happens
  //at the MFP calc point at the beginning of the next step, but it's nice to
  //have it here so we can use the new lattice info to help figure out vdir,
  //etc.) 
  this->SetLattice(postLattice);
  UpdateSCAfterLatticeChange();

  //Now, after the lattice update, we compute the post-volume rates for the
  //new lattice. For the ones that require energy dependence, we have to do
  //something a bit janky to get the updated energy dependence. If we're
  //entering a lattice we have never entered before, the rate-vs-energy array
  //has to be calculated. A few notes:
  //1. Usually this happens in GetMFP for that process, where the process has
  //   full permissions to update itself and its data members (i.e. rate models)
  //2. However, we're outside of all of those classes, accessing them via the
  //   processManager, which only allows access through const functions which
  //   don't allow us to update the "longitudinally-used" rate model for
  //   this particle. Hence, can't "just" call rate model rate.
  //3. To circumvent this problem, we accept a small CPU hit, and create a
  //   new rate model function, fill it with the new lattice, and query ITS
  //   rate. It then gets destroyed but we have the information, which is what
  //   matters. Once the subsequent turnaround (or transmission) step happens
  //   the real process will have the next lattice in its lookup table map,
  //   which means that if a QP is the first thing to enter a physical lattice,
  //   we duplicate work a bit.
  //   a.) Note that in the scenario where we move away from the "compute as
  //       needed" lookup tables, (as we have said we're going to do in a
  //       ticket) then this condition will nominally be eliminated.
  //4. However, if the Lattice info has already been calculated, then we just
  //   use the lookup table's entry for that map. Trivial.
  //Again, apologies for making this kind of janky, and I too hate this with
  //all of my soul. But this adaptive bias correction, much like the entire
  //Walk-on-spheres diffusion algorithm itself, is far enough outside of what
  //G4 is meant to do that doing this this way is the only way I can think to
  //do it...

  G4double postVolRates[nRelevProcs] = {-1,-1,-1,-1};

  // Phonon radiation
  if (relevantQPProcessesActive[0]) {

    //Check if rateModel for phonon radiation has a lookup table for this lat.
    //If not, do the janky thing
    if (dynamic_cast<G4CMPQPRadiatesPhononRate*>
        (dynamic_cast<G4CMPQPRadiatesPhononProcess*>
         (pl[relevantQPProcessIDs[0]])->GetRateModel())
        ->CheckLookupTableForLat(postLattice) ) {

      postVolRates[0] =
        dynamic_cast<G4CMPQPRadiatesPhononRate*>
        (dynamic_cast<G4CMPQPRadiatesPhononProcess*>
         (pl[relevantQPProcessIDs[0]])->GetRateModel())
        ->Rate(aTrack,postLattice);
    } else {
      
      //Create a new rate model object, have it generate its info for the
      //new lattice set above. Have to do a few operations to load data in
      //successfully, but since those functions only talk to this object,
      //the "longitudinal" rate model that follows the QP shouldn't be affected.
      G4CMPQPRadiatesPhononRate tempQPRadiatesPhononRateObj;
      if (GetCurrentTrack()){
        tempQPRadiatesPhononRateObj.LoadDataForTrack(GetCurrentTrack());
        tempQPRadiatesPhononRateObj.LoadLatticeInfoIntoSCUtils(postLattice);
        tempQPRadiatesPhononRateObj.UpdateLookupTable(postLattice);
      }
      postVolRates[0] = tempQPRadiatesPhononRateObj.Rate(aTrack);
    }
  }

  // Recombination
  if (relevantQPProcessesActive[1]) {

    //Check if rateModel for phonon radiation has a lookup table for this lat.
    //If not, do the janky thing
    if (dynamic_cast<G4CMPQPRecombinationRate*>
        (dynamic_cast<G4CMPQPRecombinationProcess*>
         (pl[relevantQPProcessIDs[1]])->GetRateModel())
        ->CheckLookupTableForLat(postLattice) ) {

      postVolRates[1] =
        dynamic_cast<G4CMPQPRecombinationRate*>
        (dynamic_cast<G4CMPQPRecombinationProcess*>
         (pl[relevantQPProcessIDs[1]])->GetRateModel())
        ->Rate(aTrack,postLattice);
    } else {
      
      //Create a new rate model object, have it generate its info for the
      //new lattice set above. Have to do a few operations to load data in
      //successfully, but since those functions only talk to this object,
      //the "longitudinal" rate model that follows the QP shouldn't be affected.
      G4CMPQPRecombinationRate tempQPRecombinationRateObj;
      if (GetCurrentTrack()){
        tempQPRecombinationRateObj.LoadDataForTrack(GetCurrentTrack());
        tempQPRecombinationRateObj.LoadLatticeInfoIntoSCUtils(postLattice);
        tempQPRecombinationRateObj.UpdateLookupTable(postLattice);
      }
      postVolRates[1] = tempQPRecombinationRateObj.Rate(aTrack);
    }
  }

  // QPDiffusionTimeStepper
  // This one is not as hard, since there is not energy dependence. Can pull
  // the value directly from the lattice information, and don't need to
  // do any muckery with dynamic casting.
  if (relevantQPProcessesActive[2]) {
    postVolRates[2] = 1.0 / postLattice->GetSCQPDiffusionStepTau();
  }

  // QPLocalTrapping
  // This one is also energy-independent (for now). Can pull the value directly
  // from the lattice information
  if (relevantQPProcessesActive[3]) {
    postVolRates[3] = 1.0 / postLattice->GetSCQPLocalTrappingTau();
  }

  //Okay, now that we have the rates, we proceed to steps 2,3 of the ovearching
  //narrative of this function (given up top), and compute average diffusion-
  //limited distances for all processes on the pre- and post-interface volumes.
  //To do this we'll need the diffusion constant for each side.
  G4double preVol_Dn = preLattice->GetSCDn();
  G4double postVol_Dn = postLattice->GetSCDn();
  G4double preVol_Erat =
    preLattice->GetSCDelta0() / aTrack.GetKineticEnergy();
  G4double postVol_Erat =
    postLattice->GetSCDelta0() / aTrack.GetKineticEnergy();
  G4double preVol_D = preVol_Dn * sqrt(1-preVol_Erat*preVol_Erat);
  G4double postVol_D = postVol_Dn * sqrt(1-postVol_Erat*postVol_Erat);

  //Now do the competition for the pre-step volume, computing displacements
  //from rates using 1/rate plugged into the diffusion equation
  double pre_disp_winner = DBL_MAX;    
  double pre_disp_geom = preStepSafety;  
  double pre_disp_qprad =
    (preVolRates[0] < 0 ? DBL_MAX : sqrt(2*preVol_D*(1.0/preVolRates[0])));  
  double pre_disp_qprec =
    (preVolRates[1] < 0 ? DBL_MAX : sqrt(2*preVol_D*(1.0/preVolRates[1])));  
  double pre_disp_qpdts =
    (preVolRates[2] < 0 ? DBL_MAX : sqrt(2*preVol_D*(1.0/preVolRates[2])));
  double pre_disp_qpltr =
    (preVolRates[3] < 0 ? DBL_MAX : sqrt(2*preVol_D*(1.0/preVolRates[3])));  
  if (pre_disp_geom < pre_disp_winner ) pre_disp_winner = pre_disp_geom;
  if (pre_disp_qprad < pre_disp_winner) pre_disp_winner = pre_disp_qprad;
  if (pre_disp_qprec < pre_disp_winner) pre_disp_winner = pre_disp_qprec;
  if (pre_disp_qpdts < pre_disp_winner) pre_disp_winner = pre_disp_qpdts;
  if (pre_disp_qpltr < pre_disp_winner) pre_disp_winner = pre_disp_qpltr;

  //Now do competition for the post-step volume, computing displacements
  //from rates using 1/rate plugged into the diffusion equation
  double post_disp_winner = DBL_MAX;    
  double post_disp_geom = postStepSafety;  
  double post_disp_qprad =
    (postVolRates[0] < 0 ? DBL_MAX : sqrt(2*postVol_D*(1.0/postVolRates[0])));  
  double post_disp_qprec =
    (postVolRates[1] < 0 ? DBL_MAX : sqrt(2*postVol_D*(1.0/postVolRates[1])));  
  double post_disp_qpdts =
    (postVolRates[2] < 0 ? DBL_MAX : sqrt(2*postVol_D*(1.0/postVolRates[2])));
  double post_disp_qpltr =
    (postVolRates[3] < 0 ? DBL_MAX : sqrt(2*postVol_D*(1.0/postVolRates[3])));  
  if (post_disp_geom < post_disp_winner ) post_disp_winner = post_disp_geom;
  if (post_disp_qprad < post_disp_winner) post_disp_winner = post_disp_qprad;
  if (post_disp_qprec < post_disp_winner) post_disp_winner = post_disp_qprec;
  if (post_disp_qpdts < post_disp_winner) post_disp_winner = post_disp_qpdts;
  if (post_disp_qpltr < post_disp_winner) post_disp_winner = post_disp_qpltr;

  //Now we have "winning" processes' diffusion circle length scales for the two
  //sides to the interface. For each of them we compute the number of returns
  //to the boundary using our heuristic MC-informed relationship.

  //This requires us to take into account a potential "microphysical" scattering
  //mean free path. For now we'll just use our Dn on each side to compute that
  //mean free path, and assume 2D transport:
  double pre_vg = 2e6 * CLHEP::m / CLHEP::s; //Hardcoded for now--> FIX!  
  double pre_lambda = 2 * preVol_Dn / pre_vg;
  double post_vg = 2e6 * CLHEP::m / CLHEP::s; //Hardcoded for now--> FIX!  
  double post_lambda = 2 * postVol_Dn / post_vg;

  //Now we use this with our winning circle lengths to compute our dimensionless
  //MFP/R_WoS "return parameter"...
  double pre_returnParam = pre_lambda / pre_disp_winner;   
  double post_returnParam = post_lambda / post_disp_winner;

  //...which in turn is used to compute the average number of attempted returns
  //to the boundary, <N>
  double pre_attReturns = 0.4 / pre_returnParam;
  double post_attReturns = 0.4 / post_returnParam;

  //Now we compute our statistical probability for reflection.
  //If pre_attemptedReturns is larger, then it means that we need to favor
  //transmission. Reflect with 0% probability, unless the user has specified
  //an additional reflection coefficient. If so, reflect with that.
  if (pre_attReturns > post_attReturns) {
    G4double reflProb = GetMaterialProperty("reflProb");
    return (G4UniformRand() < reflProb);
  }
  //Otherwise, we're coming from the smaller-<N> side, and need to favor
  //reflection. Here, we modulate user reflection by saying:
  //1. Run proper dynamic reflection, and if you DON'T reflect...
  //2. Try again, reflecting with user-specified reflection probability  
  else{
    if (G4UniformRand() < (1-pre_attReturns/post_attReturns) ) return true;
    else{
      return( G4UniformRand() < GetMaterialProperty("reflProb") );
    }
  }
  
  //For now, let's just print these and see what we get
  if (verboseLevel > 5) {
    G4cout << "In ReflectTrack: "
           << "\n pos: " << pos
           << ",\n pre-vol norm: " << pre_vol_norm
           << ",\n post-vol norm: " << post_vol_norm
           << ",\n pre-momentumDir: "
           << aStep.GetPreStepPoint()->GetMomentumDirection()
           << "\n, preStepSafety: " << preStepSafety
           << "\n, postStepSafety: " << postStepSafety
           << "\n, pre-vol qpRadiatesPhonon rate: " << preVolRates[0]
           << "\n, post-vol qpRadiatesPhonon rate: " << postVolRates[0]
           << "\n, pre-vol qpRecombination rate: " << preVolRates[1]
           << "\n, post-vol qpRecombination rate: " << postVolRates[1]
           << "\n, pre-vol qpDiffusionTimeStepper rate: " << preVolRates[2]
           << "\n, post-vol qpDiffusionTimeStepper rate: " << postVolRates[2]
           << "\n, pre-vol qpLocalTrapping rate: " << preVolRates[3]
           << "\n, post-vol qpLocalTrapping rate: " << postVolRates[3]    
           << G4endl;   
  }
}


// Do absorption of a quasiparticle. The way this is currently run is
// to pass the buck to the BoundaryUtils class, which does a NIEL
// calculation and passage to the partitioner. If we want QPs not to
// go through that process, then we can switch to just the final
// two lines of this function that are commented out.
void G4CMPQPBoundaryProcess::DoAbsorption(const G4Track& aTrack,
                                          const G4Step& aStep,
                                          G4ParticleChange& /*particleChange*/) {
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
                                          G4ParticleChange& /*particleChange*/) {

  //Debugging
  if (verboseLevel > 5) {
    G4cout << "-- G4CMPQPBoundaryProcess::DoReflection() --" << G4endl;
    G4cout << "RWBoundary DR Function Point A | Using reflection where all "
           << "returned directions are surface norms." << G4endl;
  }
  
  //This function is to be used with QP diffusion. It *will* return the
  //momentum as the surface normal in the direction of motion, as that
  //information is needed by/used by the diffusion class.
  if (verboseLevel>1) {
    G4cout << procName << ": Track reflected "
           << G4CMP::GetTrackInfo<G4CMPVTrackInfo>(aTrack)->ReflectionCount()
           << " times." << G4endl;
  }
  
  //To determine new random direction, need to know relationship between current
  //direction and the surface normal. If they are more parallel, then we need
  //to ensure that the new direction dotted into the norm is negative. If they
  //are more antiparallel, we need to make sure that the new direction dotted into the norm is positive.
  G4ThreeVector pdir = aTrack.GetMomentumDirection();
  G4ThreeVector norm = G4CMP::GetSurfaceNormal(aStep,pdir);

  //WIth the norm now a generalized surface norm, the return direction should
  //be in the opposite direction of that norm. The logic is handled
  //within GetSurfaceNormal.
  G4ThreeVector newDir = -1*norm; 
  
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
//REL: We need to be absolutely sure that nothing here depends on the old,
//i.e. pre-step lattice, because at this point the lattice info has been
//updated
void G4CMPQPBoundaryProcess::DoTransmission(const G4Track& aTrack,
                                            const G4Step& aStep,
                                            G4ParticleChange& /*particleChange*/) {
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

  //We'll note here that the lattice changeover that USED to be here
  //is now be in ReflectTrack(), which always runs unless an absorption
  //happens.  
  G4ThreeVector vdir = aTrack.GetMomentumDirection();
  G4ThreeVector norm = G4CMP::GetSurfaceNormal(aStep,vdir);
  aParticleChange.ProposeMomentumDirection(norm); 
}

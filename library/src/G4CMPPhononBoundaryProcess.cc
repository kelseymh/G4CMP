/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4PhononReflection.cc
/// \brief Implementation of the G4PhononReflection class
//
// This process handles the interaction of phonons with
// boundaries. Implementation of this class is highly
// geometry dependent.Currently, phonons are killed when
// they reach a boundary. If the other side of the
// boundary was Al, a hit is registered.
//
// $Id$
//
// 20131115  Throw exception if track's polarization state is invalid.
// 20140103  Move charge version of code here, still commented out
// 20140331  Add required process subtype code
// 20160624  Use GetTrackInfo() accessor
// 20160903  Migrate to use G4CMPBoundaryUtils for most functionality
// 20160906  Follow constness of G4CMPBoundaryUtils
// 20161114  Use new G4CMPPhononTrackInfo
// 20170829  Add detailed diagnostics to identify boundary issues
// 20170928  Replace "pol" with "mode" for phonons
// 20181010  J. Singh -- Use new G4CMPAnharmonicDecay for boundary decays
// 20181011  M. Kelsey -- Add LoadDataForTrack() to initialize decay utility.
// 20220712  M. Kelsey -- Pass process pointer to G4CMPAnharmonicDecay
// 20220905  G4CMP-310 -- Add increments of kPerp to avoid bad reflections.
// 20220910  G4CMP-299 -- Use fabs(k) in absorption test.

#include "G4CMPPhononBoundaryProcess.hh"
#include "G4CMPAnharmonicDecay.hh"
#include "G4CMPConfigManager.hh"
#include "G4LatticeManager.hh"
#include "G4CMPGeometryUtils.hh"
#include "G4CMPPhononTrackInfo.hh"
#include "G4CMPSurfaceProperty.hh"
#include "G4CMPTrackUtils.hh"
#include "G4CMPUtils.hh"
#include "G4ExceptionSeverity.hh"
#include "G4LatticePhysical.hh"
#include "G4ParallelWorldProcess.hh"
#include "G4ParticleChange.hh"
#include "G4PhysicalConstants.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VParticleChange.hh"
#include "G4VSolid.hh"
#include "Randomize.hh"
#include "G4PhononLong.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononTransSlow.hh"
#include "G4PhononScattering.hh"
#include "G4PhononPolycrystalElasticScattering.hh"
#include "G4PhononDownconversion.hh"
#include "G4CMPTrackLimiter.hh"

// Constructor and destructor

G4CMPPhononBoundaryProcess::G4CMPPhononBoundaryProcess(const G4String& aName)
  : G4VPhononProcess(aName, fPhononReflection), G4CMPBoundaryUtils(this),
    anharmonicDecay(new G4CMPAnharmonicDecay(this)) {;}

G4CMPPhononBoundaryProcess::~G4CMPPhononBoundaryProcess() {
  delete anharmonicDecay;
}


// Configure for current track including AnharmonicDecay utility

void G4CMPPhononBoundaryProcess::LoadDataForTrack(const G4Track* track) {
  G4CMPProcessUtils::LoadDataForTrack(track);
  anharmonicDecay->LoadDataForTrack(track);
}


// Compute and return step length

G4double G4CMPPhononBoundaryProcess::
PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
                                     G4double previousStepSize,
                                     G4ForceCondition* condition) {
  return GetMeanFreePath(aTrack, previousStepSize, condition);
}

G4double G4CMPPhononBoundaryProcess::GetMeanFreePath(const G4Track& aTrack,
                                             G4double /*prevStepLength*/,
                                             G4ForceCondition* condition) {
  G4cout << "REL GetMFP_G4CMPPhononBoundaryProcess" << G4endl;
  UpdateMeanFreePathForLatticeChangeover(aTrack);
  
  *condition = Forced;
  return DBL_MAX;
}


// Process action

G4VParticleChange*
G4CMPPhononBoundaryProcess::PostStepDoIt(const G4Track& aTrack,
                                         const G4Step& aStep) {
  // NOTE:  G4VProcess::SetVerboseLevel is not virtual!  Can't overlaod it
  G4CMPBoundaryUtils::SetVerboseLevel(verboseLevel);

  aParticleChange.Initialize(aTrack);
  if (!IsGoodBoundary(aStep))
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);

  if (verboseLevel>1) G4cout << GetProcessName() << "::PostStepDoIt" << G4endl;

  if (verboseLevel>2) {
    G4cout << " K direction: " << GetLocalWaveVector(aTrack).unit()
           << "\n P direction: " << aTrack.GetMomentumDirection() << G4endl;
  }

  G4cout << "REL just about to apply boundary action" << G4endl;
  
  ApplyBoundaryAction(aTrack, aStep, aParticleChange);

  ClearNumberOfInteractionLengthLeft();		// All processes should do this!
  return &aParticleChange;
}


G4bool G4CMPPhononBoundaryProcess::AbsorbTrack(const G4Track& aTrack,
                                               const G4Step& aStep) const {
  G4double absMinK = GetMaterialProperty("absMinK");
  G4ThreeVector k = G4CMP::GetTrackInfo<G4CMPPhononTrackInfo>(aTrack)->k();

  if (verboseLevel>1) {
    G4cout << GetProcessName() << "::AbsorbTrack() k " << k
	   << "\n |k_perp| " << fabs(k*G4CMP::GetSurfaceNormal(aStep))
	   << " vs. absMinK " << absMinK << G4endl;
  }

  return (G4CMPBoundaryUtils::AbsorbTrack(aTrack,aStep) &&
	  fabs(k*G4CMP::GetSurfaceNormal(aStep)) > absMinK);
}


void G4CMPPhononBoundaryProcess::
DoReflection(const G4Track& aTrack, const G4Step& aStep,
	     G4ParticleChange& particleChange) {
  auto trackInfo = G4CMP::GetTrackInfo<G4CMPPhononTrackInfo>(aTrack);

  if (verboseLevel>1) {
    G4cout << GetProcessName() << ": Track reflected "
           << trackInfo->ReflectionCount() << " times." << G4endl;
  }

  G4ThreeVector waveVector = trackInfo->k();
  G4int mode = GetPolarization(aStep.GetTrack());
  G4ThreeVector surfNorm = G4CMP::GetSurfaceNormal(aStep);

  if (verboseLevel>2) {
    G4cout << "\n Old wavevector direction " << waveVector.unit()
	   << "\n Old momentum direction   " << aTrack.GetMomentumDirection()
	   << G4endl;
  }

  // Check whether step has proper boundary-stopped geometry
  G4ThreeVector surfacePoint;
  if (!CheckStepBoundary(aStep, surfacePoint)) {
    G4cout << "REL checking step boundary failed in DoReflection" << G4endl;
    if (verboseLevel>2)
      G4cout << " Boundary point moved to " << surfacePoint << G4endl;

    particleChange.ProposePosition(surfacePoint);	// IS THIS CORRECT?!?
  }

  G4double freq = GetKineticEnergy(aTrack)/h_Planck;	// E = hf, f = E/h
  G4double specProb = surfProp->SpecularReflProb(freq);
  G4double diffuseProb = surfProp->DiffuseReflProb(freq);
  G4double downconversionProb = surfProp->AnharmonicReflProb(freq);

  // Empirical functions may lead to non normalised probabilities.
  // Normalise here.

  G4double norm = specProb + diffuseProb + downconversionProb;

  specProb /= norm;
  diffuseProb /= norm;
  downconversionProb /= norm;

  G4ThreeVector reflectedKDir;

  G4double random = G4UniformRand();

  if (verboseLevel > 2) {
    G4cout << "Surface Downconversion Probability: " << downconversionProb
	   << " random: " << random << G4endl;
  }

  G4String refltype = "";		// For use in failure message if needed

  if (random < downconversionProb) {
    if (verboseLevel > 2) G4cout << " Anharmonic Decay at boundary." << G4endl;

    /* Do Downconversion */
    anharmonicDecay->DoDecay(aTrack, aStep, particleChange);
    G4Track* sec1 = particleChange.GetSecondary(0);
    G4Track* sec2 = particleChange.GetSecondary(1);

    G4ThreeVector vec1 = GetLambertianVector(surfNorm, mode);
    G4ThreeVector vec2 = GetLambertianVector(surfNorm, mode);

    sec1->SetMomentumDirection(vec1);
    sec2->SetMomentumDirection(vec2);

    return;
  } else if (random < downconversionProb + specProb) {
    reflectedKDir = GetReflectedVector(waveVector, surfNorm, mode);
    refltype = "specular";
  } else {
    reflectedKDir = GetLambertianVector(surfNorm, mode);
    refltype = "diffuse";
  }

  G4ThreeVector vdir = theLattice->MapKtoVDir(mode, reflectedKDir);
  G4double v = theLattice->MapKtoV(mode, reflectedKDir);

  if (verboseLevel>2) {
    G4cout << "\n New wavevector direction " << reflectedKDir
	   << "\n New momentum direction   " << vdir << G4endl;
  }

  // If reflection failed, report problem and kill the track
  if (!G4CMP::PhononVelocityIsInward(theLattice,mode,reflectedKDir,surfNorm)) {
    G4Exception((GetProcessName()+"::DoReflection").c_str(), "Boundary010",
		JustWarning, ("Phonon "+refltype+" reflection failed").c_str());
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

  trackInfo->SetWaveVector(reflectedKDir);
  particleChange.ProposeVelocity(v);
  particleChange.ProposeMomentumDirection(vdir);
}


// Generate specular reflection corrected for momentum dispersion

G4ThreeVector G4CMPPhononBoundaryProcess::
GetReflectedVector(const G4ThreeVector& waveVector,
		   const G4ThreeVector& surfNorm, G4int mode) const {
  // Specular reflecton should reverses momentum along normal
  G4ThreeVector reflectedKDir = waveVector.unit();
  G4double kPerp = reflectedKDir * surfNorm;
  (reflectedKDir -= 2.*kPerp*surfNorm).setMag(1.);
  
  if (verboseLevel>2) {
    G4cout << " specular reflection with normal " << surfNorm
	   << "\n Perpendicular wavevector " << kPerp*surfNorm
	   << " (mag " << kPerp << ")" << G4endl;
  }
  
  if (G4CMP::PhononVelocityIsInward(theLattice,mode,reflectedKDir,surfNorm))
    return reflectedKDir;

  // Reflection didn't work as expected, need to correct   

  // Watch how momentum direction changes with each kPerp step
  G4ThreeVector olddir, newdir;
  
  olddir = theLattice->MapKtoVDir(mode, reflectedKDir);
  G4double kstep = 0.1*kPerp;
  G4int nstep = 0, nflip = 0;
  while (fabs(kstep) > 1e-6 && fabs(nstep*kstep)<1. && 
	 !G4CMP::PhononVelocityIsInward(theLattice,mode,reflectedKDir,surfNorm)) {
    newdir = theLattice->MapKtoVDir(mode, reflectedKDir);
    if (newdir*surfNorm > olddir*surfNorm && nflip<5) {
      if (verboseLevel>2) {
	G4cout << " Reflected wv pushing momentum outward:"
	       << " newdir*surfNorm = " << newdir*surfNorm
	       << G4endl;
      }
      
      kstep = -kstep;
      nflip++;
    }
    
    (reflectedKDir -= kstep*surfNorm).setMag(1.);
    olddir = newdir;
    nstep++;
  } 
  
  if (nstep>0 && verboseLevel) {
    G4cout << " adjusted specular reflection with " << nstep << " steps"
	   << " (" << nflip << " flips) kPerp " << kPerp << G4endl;
  }

  return reflectedKDir;
}


// Generate diffuse reflection according to 1/cos distribution

G4ThreeVector G4CMPPhononBoundaryProcess::
GetLambertianVector(const G4ThreeVector& surfNorm, G4int mode) const {
  G4ThreeVector reflectedKDir;
  const G4int maxTries = 1000;
  G4int nTries = 0;
  do {
    reflectedKDir = G4CMP::LambertReflection(surfNorm);
  } while (nTries++ < maxTries &&
	   !G4CMP::PhononVelocityIsInward(theLattice, mode,
					  reflectedKDir, surfNorm));

  return reflectedKDir;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//Now that we can actually mate two different materials and have phonons propagate
//between them, we should make a nontrivial phonon doTransmission class.
void G4CMPPhononBoundaryProcess::DoTransmission(const G4Track& aTrack,
						const G4Step& aStep,
						G4ParticleChange& aParticleChange) {

  //Generic print
  if (verboseLevel>1) G4cout << "REL-- Track transmission requested" << G4endl;

  G4cout << "DoTransmision: aTrack getposition: " << aTrack.GetPosition() << G4endl;
  G4cout << "aStep poststepposition: " << aStep.GetPostStepPoint()->GetPosition() << G4endl;
  G4cout << "Lattice manager current lattice: " << G4LatticeManager::GetLatticeManager()->GetLattice(aStep.GetPreStepPoint()->GetPhysicalVolume()) << ", lattice manager post step point lattice: " << G4LatticeManager::GetLatticeManager()->GetLattice(aStep.GetPostStepPoint()->GetPhysicalVolume()) << G4endl;
  
  //First, check the different materials involved. If one of them does not have a lattice, then something upstream has gone wrong. Scream.
  //Note that here, we need to use the G4LatticeManager rather than the TrackUtils version of GetLattice() because the trackUtils
  //version somehow chooses a point that's inside the current (i.e. not the post-step) volume.
  G4LatticePhysical * latNear = G4LatticeManager::GetLatticeManager()->GetLattice(aStep.GetPreStepPoint()->GetPhysicalVolume());
  G4LatticePhysical * latFar = G4LatticeManager::GetLatticeManager()->GetLattice(aStep.GetPostStepPoint()->GetPhysicalVolume());
  if( !latNear || !latFar ){
    G4ExceptionDescription msg;
    msg << "Expecting to do phonon transmission at interface but one or more lattices splitting the boundary cannot be found.";
    G4Exception("G4CMPPhononBoundaryProcess::DoTransmission", "PhononBoundary001",FatalException,msg);
  }

  //SIMPLEST POSSIBLE IMPLEMENTATION -- PHYSICS IS NOT NECESSARILY RIGHT
  //Now get track info at this point: the k vector, mode, etc.
  auto trackInfo = G4CMP::GetTrackInfo<G4CMPPhononTrackInfo>(aTrack);
  G4ThreeVector waveVector = trackInfo->k();
  G4int mode = GetPolarization(aStep.GetTrack());
  
  // Check whether step has proper boundary-stopped geometry
  G4ThreeVector surfacePoint;
  if (!CheckStepBoundary(aStep, surfacePoint)) {
    G4cout << "REL checking step boundary failed in DoTransmission" << G4endl;
    if (verboseLevel>2)
      G4cout << " Boundary point moved to " << surfacePoint << G4endl;

    aParticleChange.ProposePosition(surfacePoint);	// IS THIS CORRECT?!?
  }

  
  //Since the lattice hasn't changed yet, change it here. (This also happens at the MFP calc point at the beginning of the next step,
  //but it's nice to have it here so we can use the new lattice info to help figure out vdir, etc.)
  this->SetLattice(G4LatticeManager::GetLatticeManager()->GetLattice(aStep.GetPostStepPoint()->GetPhysicalVolume()));
  UpdateSCAfterLatticeChange();

  G4cout << "REL the lattice at the end of doTransmission: " << theLattice << G4endl;
  G4ThreeVector surfNorm = G4CMP::GetSurfaceNormal(aStep);
  G4ThreeVector vdir = theLattice->MapKtoVDir(mode, waveVector);
  G4double v = theLattice->MapKtoV(mode, waveVector);
  trackInfo->SetWaveVector(waveVector);
  aParticleChange.ProposeVelocity(v);
  aParticleChange.ProposeMomentumDirection(vdir);

  
}

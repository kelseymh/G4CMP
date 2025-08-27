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
// 20240718  G4CMP-317 -- Initial implementation of surface displacement.
// 20250124  G4CMP-447 -- Use FillParticleChange() to update wavevector and Vg.
// 20250204  G4CMP-459 -- Support reflection displacement search at hard corners.
// 20250325  G4CMP-463 -- Set surface step size & limit with macro command.
// 20250402  G4CMP-468 -- Support position change after surface displacement.
// 20250413  M. Kelsey -- Protect debugging output with verbosity.
// 20250422  N. Tenpas -- Add position arguments for PhononVelocityIsInward test.
// 20250423  G4CMP-468 -- Add wrapper function for updating navigator.
// 20250423  G4CMP-468 -- Move GetLambertianVector to G4CMPUtils.
// 20250424  G4CMP-465 -- Move custom solid functions to new G4CMPSolidUtils.
// 20250429  G4CMP-461 -- Implement ability to skip flats during displacement.
// 20250505  G4CMP-458 -- Rename GetReflectedVector to GetSpecularVector.
// 20250505  G4CMP-471 -- Update diagnostic output for surface displacement loop.
// 20250814  G4CMP-496 -- Fix phonon wavevector & momentum direction coming out of DoDecay.

#include "G4CMPPhononBoundaryProcess.hh"
#include "G4CMPAnharmonicDecay.hh"
#include "G4CMPConfigManager.hh"
#include "G4LatticeManager.hh"
#include "G4CMPGeometryUtils.hh"
#include "G4CMPParticleChangeForPhonon.hh"
#include "G4CMPPhononTrackInfo.hh"
#include "G4CMPSolidUtils.hh"
#include "G4CMPSurfaceProperty.hh"
#include "G4CMPTrackUtils.hh"
#include "G4CMPUtils.hh"
#include "G4ExceptionSeverity.hh"
#include "G4GeometryTolerance.hh"
#include "G4LatticePhysical.hh"
#include "G4Navigator.hh"
#include "G4ParallelWorldProcess.hh"
#include "G4ParticleChange.hh"
#include "G4PhononPolarization.hh"
#include "G4PhysicalConstants.hh"
#include "G4RandomDirection.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4TransportationManager.hh"
#include "G4UnitsTable.hh"
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
#include <float.h>

// Constructor and destructor

G4CMPPhononBoundaryProcess::G4CMPPhononBoundaryProcess(const G4String& aName)
  : G4VPhononProcess(aName, fPhononReflection), G4CMPBoundaryUtils(this),
    anharmonicDecay(new G4CMPAnharmonicDecay(this)), stepSize(0*um), nStepLimit(0) {
  // Initialize stepSize and max step limit from config manager
  G4CMPConfigManager* config = G4CMPConfigManager::Instance();
  stepSize = config->GetPhononSurfStepSize();
  nStepLimit = config->GetPhononSurfStepLimit();
  // Register custom ParticleChange with G4VProcess base class
  pParticleChange = &phParticleChange;
}

G4CMPPhononBoundaryProcess::~G4CMPPhononBoundaryProcess() {
  delete anharmonicDecay;
}


// Configure for current track including AnharmonicDecay utility

void G4CMPPhononBoundaryProcess::LoadDataForTrack(const G4Track* track, const G4bool overrideMomentumReset) {
  //G4cout << "In G4CMPPhononBoundaryProcess::LoadDataForTrack(), "
  //<< " running the ProcessUtils version." << G4endl;
  G4CMPProcessUtils::LoadDataForTrack(track,overrideMomentumReset);

  //G4cout << "In G4CMPPhononBoundaryProcess::LoadDataForTrack(), "
  //<< " running the Anharmonic Decay version." << G4endl;
  anharmonicDecay->LoadDataForTrack(track,overrideMomentumReset);
}


// Compute and return step length

G4double G4CMPPhononBoundaryProcess::
PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
                                     G4double previousStepSize,
                                     G4ForceCondition* condition) {
  return GetMeanFreePath(aTrack, previousStepSize, condition);

  if( verboseLevel > 5 ){
    G4cout << "In PostStepGetPhysicalInteractionLength: Track momentum"
	   << "direction: " << aTrack.GetMomentumDirection() << G4endl;
    G4cout << "In PostStepGetPhysicalInteractionLength: Track position: "
	   << aTrack.GetPosition() << G4endl;
  }  
}

G4double G4CMPPhononBoundaryProcess::GetMeanFreePath(const G4Track& aTrack,
                                             G4double /*prevStepLength*/,
                                             G4ForceCondition* condition) {
  if (verboseLevel > 5) {
    G4cout << "-- G4CMPPhononBoundaryProcess::GetMeanFreePath --" << G4endl;
    G4cout << "Track momentum direction: " << aTrack.GetMomentumDirection()
	   << G4endl;
  }
  UpdateMeanFreePathForLatticeChangeover(aTrack);

  //G4cout << "REL Lattice seen by anharmonic decay is "
  //<< anharmonicDecay->GetLattice()->GetLattice()->GetName() << G4endl;
  
  if (verboseLevel > 5) {
    G4cout << "Post UpdateMeanFreePathForLatticeChamgeover: Track momentum "
	   << "direction: " << aTrack.GetMomentumDirection() << G4endl;
    G4cout << "Post UpdateMeanFreePathForLatticeChangeover: Track position: "
	   << aTrack.GetPosition() << G4endl;
  }
  
  *condition = Forced;
  return DBL_MAX;
}


// Process action

G4VParticleChange*
G4CMPPhononBoundaryProcess::PostStepDoIt(const G4Track& aTrack,
                                         const G4Step& aStep) {
  // NOTE:  G4VProcess::SetVerboseLevel is not virtual!  Can't overlaod it
  G4CMPBoundaryUtils::SetVerboseLevel(verboseLevel);  
  
  //Debugging
  if( verboseLevel > 5 ){
    G4cout << "-- G4CMPPhononBoundaryProcess::PostStepDoIt --" << G4endl;
    G4cout << "Track momentum direction: " << aTrack.GetMomentumDirection()
	   << G4endl;
    G4cout << "Track position: " << aTrack.GetPosition() << G4endl;
    G4cout << "Step pre-step point: " << aStep.GetPreStepPoint()->GetPosition()
	   << G4endl;
    G4cout << "Step post-step point: "
	   << aStep.GetPostStepPoint()->GetPosition() << G4endl;
    G4cout << "Step pre-step momentum direction: "
	   << aStep.GetPreStepPoint()->GetMomentumDirection() << G4endl;
    G4cout << "Step post-step momentum direction: "
	   << aStep.GetPostStepPoint()->GetMomentumDirection() << G4endl;
  }  
  
  phParticleChange.Initialize(aTrack);
  

  if (!IsGoodBoundary(aStep))
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);

  if (verboseLevel>1) {
    G4int eID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    G4cout << GetProcessName() << "::PostStepDoIt "
           << "Event " << eID << " Track " << aTrack.GetTrackID()
	   << " Step " << aTrack.GetCurrentStepNumber() << G4endl;
  }

  if (verboseLevel>2) {
    G4cout << " K direction: " << GetLocalWaveVector(aTrack).unit()
           << "\n P direction: " << aTrack.GetMomentumDirection() << G4endl;
  }


  //Debugging
  if( verboseLevel > 5 ){
    G4cout << "PSDI Function Point A | just about to apply boundary action"
	   << G4endl;
  }
  
  ApplyBoundaryAction(aTrack, aStep, phParticleChange);
  //ClearNumberOfInteractionLengthLeft(); // All processes should do this!
  //REL This^ was removed by NT or WL, I think? Why?
  
  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
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

  //Debugging
  if( verboseLevel > 5 ){
    G4cout << "-- G4CMPPhononBoundaryProcess::DoReflection --" << G4endl;
    G4cout << "DR Function Point A | aStep pre-step touchable: "
	   << aStep.GetPreStepPoint()->GetTouchable()->GetVolume()->GetName()
	   << G4endl;
    G4cout << "DR Function Point A | aStep post-step touchable: "
	   << aStep.GetPostStepPoint()->GetTouchable()->GetVolume()->GetName()
	   << G4endl;
    G4cout << "DR Function Point A | current touchable: "
	   << GetCurrentTouchable()->GetVolume()->GetName() << G4endl;
  }

  
  auto trackInfo = G4CMP::GetTrackInfo<G4CMPPhononTrackInfo>(aTrack);

  if (verboseLevel>1) {
    G4cout << GetProcessName() << ": Track reflected "
           << trackInfo->ReflectionCount() << " times." << G4endl;
  }

  G4ThreeVector waveVector = trackInfo->k();
  G4ThreeVector initialVDir = aStep.GetPreStepPoint()->GetMomentumDirection();
  if (verboseLevel > 5) {
    G4cout << "DR Function Point B | initialVDir check: " << initialVDir
	   << G4endl;
  }
  
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

    //Debugging
    if (verboseLevel > 5) {
      G4cout << "DR Function Point C | checking step boundary failed in "
	     << "DoReflection" << G4endl;
    }
 
    if (verboseLevel>2)
      G4cout << " Boundary point moved to " << surfacePoint << G4endl;

    particleChange.ProposePosition(surfacePoint);	// IS THIS CORRECT?!?
  }

  G4double freq = GetKineticEnergy(aTrack)/h_Planck;	// E = hf, f = E/h
  G4double specProb = surfProp->SpecularReflProb(freq);
  G4double diffuseProb = surfProp->DiffuseReflProb(freq);
  G4double downconversionProb = surfProp->AnharmonicReflProb(freq);
  //G4cout << "AnharmonicProb at point: " << surfacePoint << " is "
  //	 << downconversionProb << G4endl;
  //G4cout << "DiffuseProb at point: " << surfacePoint << " is "
  //	 << diffuseProb << ", spec: " << specProb << G4endl;
  //G4cout << "SpecProb at point: " << surfacePoint << " is "
  //	 << specProb << G4endl;
  
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

    //Debugging
    if( verboseLevel > 5 ){
      G4cout << "DR Function Point CA | after DoDecay in surface-mediated "
	     << "downconversion" << G4endl;
    }
  
    G4int mode1 = G4PhononPolarization::Get(sec1->GetParticleDefinition());
    G4int mode2 = G4PhononPolarization::Get(sec2->GetParticleDefinition());

    G4ThreeVector vec1 = G4CMP::GetLambertianVector(theLattice, surfNorm, mode1,
                                                    surfacePoint);
    G4ThreeVector vec2 = G4CMP::GetLambertianVector(theLattice, surfNorm, mode2,
                                                    surfacePoint);

    UpdatePhononWavevector(*sec1, vec1);
    UpdatePhononWavevector(*sec2, vec2);

    //Debugging
    if( verboseLevel > 5 ){
      G4cout << "DR Function Point CB | after GetPolarizations in "
	     << "surface-mediated downconversion" << G4endl;
    }

    //Keep in mind that the surfNorm here is in global frame. Here, the input
    //vDir to be that of the incident, singular phonon because all we need it
    //for is orientation.
    G4ThreeVector vec1 = G4CMP::GetLambertianVector(theLattice, surfNorm, mode1,
						    initialVDir, surfacePoint);
    G4ThreeVector vec2 = G4CMP::GetLambertianVector(theLattice, surfNorm, mode2,
						    initialVDir,surfacePoint);

    //Catch exceptions
    if( vec1.mag() == 0 || vec2.mag() == 0 ){
      G4String msg = "Either vec1 or vec2 in downconversion-at-a-boundary has "
	"zero length, implying that one of the Lambertian vector attempts has "
	"failed. Here we're going to overdo it and kill the whole program, but "
	"should come back later to figure out a more efficient way to kill "
	"these tracks.";
      G4Exception((GetProcessName()+"::DoReflection").c_str(), "Boundary012",
		  FatalException, msg.c_str());
      return;
    }

      
    
    //Debugging
    if (verboseLevel > 5) {
      G4cout << "DR Function Point D | In anharmonic decay, wavevector 1 just "
	     << "post-calculate is: " << vec1 << G4endl;
      G4cout << "DR Function Point D | In anharmonic decay, wavevector 2 just "
	     << "post-calculate is: " << vec2 << G4endl;
    }
    
    //We've now got lambertian reflected k-vectors in the global frame. Need
    //to do two things:
    //1. Store this k-vector info in the track info
    auto trackInfo1 = G4CMP::GetTrackInfo<G4CMPPhononTrackInfo>(sec1);
    auto trackInfo2 = G4CMP::GetTrackInfo<G4CMPPhononTrackInfo>(sec2);
    trackInfo1->SetWaveVector(vec1);
    trackInfo2->SetWaveVector(vec2);
    
    //2. Generate a direction from this. This requires converting the global
    //waveVector to a local direction, converting to velocity, and then
    //converting velocity back to global direction
    RotateToLocalDirection(vec1);
    RotateToLocalDirection(vec2);

    //Debugging
    if (verboseLevel > 5){
      G4cout << "DR Function Point E | In anharmonic decay, wavevector 1 "
	     << "rotated to local direction is: " << vec1 << G4endl;
      G4cout << "DR Function Point E | In anharmonic decay, wavevector 2 "
	     << "rotated to local direction is: " << vec2 << G4endl;
    }
    G4ThreeVector vDir1 = theLattice->MapKtoVDir(mode1, vec1);    
    G4ThreeVector vDir2 = theLattice->MapKtoVDir(mode2, vec2);

    //Debugging
    if (verboseLevel > 5){
      G4cout << "DR Function Point F | In anharmonic decay, vDir 1 in local "
	     << "direction is: " << vDir1 << G4endl;
      G4cout << "DR Function Point F | In anharmonic decay, vDir 2 in local "
	     << "direction is: " << vDir2 << G4endl;
    }
    
    RotateToGlobalDirection(vDir1);
    RotateToGlobalDirection(vDir2);
    sec1->SetMomentumDirection(vDir1);
    sec2->SetMomentumDirection(vDir2);

    //Debugging
    if (verboseLevel > 5){      
      G4cout << "DR Function Point G | In anharmonic decay, wavevector 1 just "
	     << "pre-return is: " << vec1 << " corresponds to global velocity "
	     << "direction: " << vDir1 << G4endl;      
      G4cout << "DR Function Point G | In anharmonic decay, wavevector 2 just "
	     << "pre-return is: " << vec2 << " corresponds to global velocity "
	     << "direction: " << vDir2 << G4endl;
    }
    
    return;
  } else if (random < downconversionProb + specProb) {

    if (verboseLevel > 2) {
      G4cout << "Specular reflection at boundary, over surface norm "
	     << surfNorm << "." << G4endl;
    }
    
    // Modify surfacePoint & surfNorm in place
    reflectedKDir = GetSpecularVector(waveVector, surfNorm, mode, initialVDir,
				      surfacePoint);
    
    refltype = "specular";
  } else {
    if (verboseLevel > 2 ){
      G4cout << "Diffuse reflection at boundary, with surface norm "
	     << surfNorm << ", in lattice: " << theLattice
	     << ", at surfacePoint: " << surfacePoint << "." << G4endl;
    }

    //Debugging
    if (verboseLevel > 5) {
      G4cout << "DR Function Point H | Diffuse reflection at boundary, with "
	     << "surface norm " << surfNorm << ", in lattice: "
	     << theLattice << ", at surfacePoint: " << surfacePoint << "."
	     << G4endl;
    }

    reflectedKDir = G4CMP::GetLambertianVector(theLattice, surfNorm, mode,
					       initialVDir, surfacePoint);    
    refltype = "diffuse";
  }

  //Do a check to make sure that we catch exceptions (here not exactly zero
  //to catch floating point exceptions)
  if( reflectedKDir.mag() < 1e-10 ){
    G4String msg = "ReflectedKDir has zero length, implying that a Lambertian "
      "vector attempt has failed. Will kill this phonon.";
    G4Exception((GetProcessName()+"::DoReflection").c_str(), "Boundary013",
		JustWarning, msg.c_str());
    DoSimpleKill(aTrack, aStep, particleChange);
  }

  // Update trackInfo wavevector and particleChange's group velocity and
  // momentum direction
  // reflectedKDir is in global coordinates here - no conversion needed
  FillParticleChange(particleChange, aTrack, reflectedKDir);
  const G4ThreeVector vdir = *particleChange.GetMomentumDirection();

  // If displacement occured: update the particle change and navigator for
  // volume assignment
  if (refltype == "specular" && *particleChange.GetPosition() != surfacePoint) {
    FillParticleChange(phParticleChange, aStep, surfacePoint);
    UpdateNavigatorVolume(aStep, surfacePoint, vdir);
  }

  //Debugging
  if (verboseLevel>2) {
    G4cout << "New surface position " << *particleChange.GetPosition()/mm 
	   << " mm\n New wavevector direction " << reflectedKDir
	   << "\n New momentum direction " << vdir << G4endl;
  }

  // If reflection failed, report problem and kill the track.
  // For this check, need to make sure that surfNorm points opposite to the
  // direction of the incident k-Vector
  if (verboseLevel > 5){
    G4cout << "DR Function Point I | Phonon incident waveVector: "
	   << waveVector << ", non-generalized surfNorm: " << surfNorm
	   << ", and attempted reflected vector " << reflectedKDir << G4endl;
  }

  //REL changed to <= 0.0 8/25/25
  if (surfNorm.dot(initialVDir.unit()) <= 0.0 ) { surfNorm *= -1; }
  if (!G4CMP::PhononVelocityIsInward(theLattice,mode,reflectedKDir,surfNorm,
				     surfacePoint)) {
    G4String msg = G4PhononPolarization::Name(mode) + " " + refltype
      + " reflection failed";

    //Debugging
    if (verboseLevel > 5) {
      G4cout << "DR Function Point J | In final fail: Phonon incident "
	     << "waveVector: " << waveVector.unit()
	     << ", generalized surfNorm: " << surfNorm
	     << ", and attempted reflected vector " << reflectedKDir << G4endl;
    }
    
    G4Exception((GetProcessName()+"::DoReflection").c_str(), "Boundary010",
		JustWarning, msg.c_str());
    DoSimpleKill(aTrack, aStep, particleChange);
    return;
  }
}

// Generate specular reflection corrected for momentum dispersion

G4ThreeVector G4CMPPhononBoundaryProcess::
GetSpecularVector(const G4ThreeVector& waveVector,
                  G4ThreeVector& surfNorm, G4int mode,
		  G4ThreeVector& initialVDir,
                  G4ThreeVector& surfacePoint) {

  //First, use our surface normal to find the generalized surface normal
  G4ThreeVector generalizedSurfNorm =
    G4CMP::GetGeneralizedSurfaceNormal(surfNorm,initialVDir);


  //Specular reflection should reverse momentum along normal. Now
  //kPerp should almost always be positive (where the "almost" is because
  //the generalized surface norm is built on relationships between the
  //surfNorm and *velocity,* not k-vector. If it is negative (because
  //we're at a glancing enough angle that the k-vector is somehow pointing
  //inward (which we have seen before), then the reflected vector should
  //inherit a "toward-surface" component, which is what would have happened
  //in the old code as well.
  G4ThreeVector reflectedKDir = waveVector.unit();
  G4double kPerp = reflectedKDir * generalizedSurfNorm;
  (reflectedKDir -= 2.*kPerp*generalizedSurfNorm).setMag(1.);
  //REL^ changed to -= 8/25/25

  //Old version
  // Specular reflecton should reverse momentum along normal
  //G4ThreeVector reflectedKDir = waveVector.unit();
  //G4double kPerp = reflectedKDir * surfNorm; // Dot product between k and norm
  //(reflectedKDir -= 2.*kPerp*surfNorm).setMag(1.); // Reflect against normal

    
  if (G4CMP::PhononVelocityIsInward(theLattice,mode,
				    reflectedKDir,generalizedSurfNorm,
                                    surfacePoint))
    return reflectedKDir;

  
  // REL: Since the "base" version of this assumes outward-pointing surface
  // normals and since our initial surfNorm should in principle point
  // outward for most SCDMS-specific volumes, we just leave the rest
  // of the code as-is and hope all is well.   
  //  Below this line, I have not changed the algorithm except for the
  // lambertian call far below, and plan on leaving this for some combination
  // of WL, NT, MK, and myself to figure out/fix, as I'm not sure how
  // generalized this is for arbitrary volumes.
  
  // Reflection didn't work as expected, need to correct:
  // If the reflected wave vector cannot propagate in the bulk
  // (i.e., the reflected k has an associated vg which is not inwardly
  // directed.)
  
  // That surface wave will propagate until it reaches a point
  // where the wave vector has an inwardly directed v

  //New version (REL)
  //RotateToLocalDirection(reflectedKDir);
  //G4ThreeVector newNorm = generalizedSurfNorm;
  //RotateToLocalDirection(newGeneralizedNorm);
  
  //Old version
  RotateToLocalDirection(reflectedKDir);
  G4ThreeVector newNorm = surfNorm;
  RotateToLocalDirection(newNorm);

  // Initialize solid object and utilities
  G4VSolid* solid = GetCurrentVolume()->GetLogicalVolume()->GetSolid();
  G4CMPSolidUtils solidUtils(solid, verboseLevel, GetProcessName());
  
  G4ThreeVector stepLocalPos = GetLocalPosition(surfacePoint);
  G4ThreeVector oldNorm = newNorm;
  G4ThreeVector oldstepLocalPos = stepLocalPos;

  ////New version (REL)
  //// Break up wavevector to perp and tan components
  //G4double kPerpMag = reflectedKDir.dot(newNorm);
  //G4ThreeVector kPerpV = kPerpMag * newNorm;		// Positive implied in kPerpMag for inward pointing (due to generalized surfNorm)
  //G4ThreeVector kTan = reflectedKDir - kPerpV;		// Get kTan: reflectedKDir = kPerpV + kTan. Think this is the same regardless because both kPerpMag and newNorm inherit minuses

  //Old version
  // Break up wavevector to perp and tan components
  G4double kPerpMag = reflectedKDir.dot(newNorm);
  G4ThreeVector kPerpV = kPerpMag * newNorm;		// Negative implied in kPerpMag for inward pointing
  G4ThreeVector kTan = reflectedKDir - kPerpV;		// Get kTan: reflectedKDir = kPerpV + kTan

  // Get axis and phi for tangent rotations
  G4ThreeVector axis = kPerpV.cross(kTan).unit();
  G4double phi = 0.;
  EInside isIn = solid->Inside(stepLocalPos);

  G4int nAttempts = 0;

  // Initialize stepSize for _this_ solid object
  G4CMPConfigManager* config = G4CMPConfigManager::Instance();
  stepSize = config->GetPhononSurfStepSize();

  // Get flat skip step size dependent on solid
  G4ThreeVector pmin(0,0,0);
  G4ThreeVector pmax(0,0,0);
  solid->BoundingLimits(pmin, pmax);
  G4double flatStepSize = (pmax - pmin).mag();

  // Set default stepSize based on solid bounding limits
  if (stepSize == 0) stepSize = flatStepSize / 1000;

  if (verboseLevel>3) {
    G4cout << "GetSpecularVector:beforeLoop -> "
	   << ", stepLocalPos = " << stepLocalPos/mm << " mm"
	   << ", reflectedKDir = " << reflectedKDir
	   << ", newNorm = " << newNorm
	   << ", kPerpMag (newNorm dot reflectedKDir) = " << kPerpMag
	   << ", kPerpV (kPerpMag * newNorm) = " << kPerpV
	   << ", kTan (reflectedKDir - kPerpV) = " << kTan
	   << ", surfaceStepSize = " << G4BestUnit(stepSize, "Length")
	   << ", nStepLimit = " << nStepLimit << G4endl;
  }

  // Assumes everything is in Global. Just add the GetGlobal in the loop conditions.
  while (!G4CMP::PhononVelocityIsInward(theLattice, mode,
   GetGlobalDirection(reflectedKDir), GetGlobalDirection(newNorm),
   GetGlobalPosition(stepLocalPos)) && nAttempts++ < nStepLimit) {
    // Save previous loop values
    oldstepLocalPos = stepLocalPos;
    oldNorm = newNorm;

    // Step along kTan direction - this point is now outside the detector
    stepLocalPos += stepSize * kTan.unit();

    //Old version
    // Get the local normal at the new surface point
    newNorm = solid->SurfaceNormal(stepLocalPos);
    // Check position status for flat skipper
    isIn = solid->Inside(stepLocalPos);

    //New version: should use a generalized surface norm -- here I'm not sure
    //what to use as the incident momentum direction here... If the
    //displacement is along a small angle, maybe we can just ask for the
    //surfNorm that, when dotted into the original generalizedSurfNorm,
    //gives a positive number? 
    
    // Get the local normal at the new surface point
    newNorm = solid->SurfaceNormal(stepLocalPos);
    // Check position status for flat skipper
    isIn = solid->Inside(stepLocalPos);

    
    // Check if the phonon is on a flat. Must be on the solid surface
    if (oldNorm == newNorm && isIn == kSurface) {
      // Adjust stepLocalPos to edge of the flat (still on the flat)
      // Modifies stepLocalPos and kTan in place
      solidUtils.AdjustOffFlats(stepLocalPos, kTan, flatStepSize, newNorm, 0);
      // Do a diffuse reflection if stuck in regression
      if (solid->Inside(stepLocalPos) != kSurface) {
        reflectedKDir = newNorm;
        break;
      }
      // Step off the flat and adjust newNorm
      stepLocalPos += stepSize * kTan.unit();
      newNorm = solid->SurfaceNormal(stepLocalPos);
    }

    // Adjust stepLocalPos back to surface of detector
    solidUtils.AdjustToClosestSurfacePoint(stepLocalPos, -newNorm);
    // Check position status for edge reflections
    isIn = solid->Inside(stepLocalPos);

    // Large normal changes and not being on surface after initial adjustment
    // indicates we are approaching an edge
    if (isIn != kSurface || newNorm * oldNorm <= 0) {
      // Reset stepLocalPos and newNorm to last valid surface point
      stepLocalPos = oldstepLocalPos;
      newNorm = oldNorm;
      // Modify stepLocalPos in place to edge position
      solidUtils.AdjustToEdgePosition(kTan, stepLocalPos, stepSize, 1);
      // Do a diffuse reflection if the adjustment failed
      if (solid->Inside(stepLocalPos) != kSurface) {
        reflectedKDir = newNorm;
        break;
      }
      // Reflect kTan against the edge - rotates & modifies kTan; modifies
      // newNorm
      solidUtils.ReflectAgainstEdge(kTan, stepLocalPos, newNorm);
    } else {
      // Rotate kTan to new position
      axis = newNorm.cross(kTan).unit();
      phi = oldNorm.azimAngle(newNorm, axis);
      kTan = kTan.rotate(axis, phi);
    }

    // Get perpendicular component of reflected k w/ new norm
    // (negative implied in kPerpMag for inward pointing)
    //REL now positive-implied?
    kPerpV = kPerpMag * newNorm;

    // Calculate new reflectedKDir (kTan + kPerpV) and Vg
    reflectedKDir = kTan + kPerpV;
    G4ThreeVector vDir = theLattice->MapKtoVDir(mode, reflectedKDir);

    if (verboseLevel>3) {
      G4cout << " GetSpecularVector:insideLoop -> "
	     << " attempts = " << nAttempts
	     << ", oldstepLocalPos = " << oldstepLocalPos/mm << " mm"
	     << ", stepLocalPos = " << stepLocalPos/mm << " mm"
	     << ", oldNorm = " << oldNorm
	     << ", newNorm = " << newNorm
	     << ", kPerpMag = " << kPerpMag
	     << ", kPerpV (kPerpMag * newNorm) = " << kPerpV
	     << ", kTan = " << kTan
	     << ", reflectedKDir (kTan + kPerpV) = " << reflectedKDir
	     << ", Phonon mode = " << G4PhononPolarization::Label(mode)
	     << ", New group velocity: " << vDir << G4endl;
    }
  }

  // Restore global coordinates to new vectors
  RotateToGlobalDirection(reflectedKDir);
  RotateToGlobalDirection(newNorm);
  RotateToGlobalPosition(stepLocalPos);

  if (!G4CMP::PhononVelocityIsInward(theLattice, mode, reflectedKDir, newNorm,
                                     stepLocalPos)) {
    if (verboseLevel) {
      G4cerr << GetProcessName() << "::GetSpecularVector"
	     << ": Phonon displacement failed after " << nAttempts - 1
	     << " attempts." << G4endl;
      if (verboseLevel>1) {
	G4cout << "Doing diffuse reflection at surface point " 
	       << surfacePoint/mm << " mm" << G4endl;
      }
    }

    // Get reflectedKDir from initial point and restore original values
    stepLocalPos = surfacePoint;
    newNorm = surfNorm;
    reflectedKDir = G4CMP::GetLambertianVector(theLattice, surfNorm, mode,
					       initialVDir, 
                                               surfacePoint);
  }

  if (verboseLevel>3) {
    G4cout << GetProcessName() << "::GetSpecularVector"
	   << ": nAttempts = " << nAttempts
	   << ", waveVector = " << waveVector
	   << ", reflectedKDir = " << reflectedKDir
	   << ", initialGlobalPostion = " << surfacePoint/mm << " mm"
	   << ", finalGlobalPosition = " << stepLocalPos/mm << " mm"
	   << G4endl;
  }

  surfacePoint = stepLocalPos;
  surfNorm = newNorm;
  return reflectedKDir;
}


void G4CMPPhononBoundaryProcess::
UpdateNavigatorVolume(const G4Step& step, const G4ThreeVector& position,
                      const G4ThreeVector& vDir) const {
  G4Navigator* navigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  navigator->LocateGlobalPointWithinVolume(position);
  G4double safety = step.GetPostStepPoint()->GetSafety();
  navigator->ComputeStep(position, vDir, step.GetStepLength(), safety);
}



//Now that we can actually mate two different materials and have phonons
//propagate between them, we should make a nontrivial phonon doTransmission
//function
void G4CMPPhononBoundaryProcess::
DoTransmission(const G4Track& aTrack,const G4Step& aStep,
	       G4ParticleChange& aParticleChange) {

  //Debugging
  if (verboseLevel > 5) {
    G4cout << "-- G4CMPPhononBoundaryProcess::DoTransmission --" << G4endl;
    G4cout << "DT Function Point A | aTrack getposition: "
	   << aTrack.GetPosition() << G4endl;
    G4cout << "DT Function Point A | aStep poststepposition: "
	   << aStep.GetPostStepPoint()->GetPosition() << G4endl;
    G4cout << "DT Function Point A | aStep pre-step touchable: "
	   << aStep.GetPreStepPoint()->GetTouchable()->GetVolume()->GetName()
	   << G4endl;
    G4cout << "DT Function Point A | aStep post-step touchable: "
	   << aStep.GetPostStepPoint()->GetTouchable()->GetVolume()->GetName()
	   << G4endl;
    G4cout << "DT Function Point A | current touchable: "
	   << GetCurrentTouchable()->GetVolume()->GetName() << G4endl;
    G4cout << "DT Function Point A | Lattice manager current lattice: "
	   << G4LatticeManager::GetLatticeManager()->GetLattice(aStep.GetPreStepPoint()->GetPhysicalVolume())
	   << ", lattice manager post step point lattice: "
	   << G4LatticeManager::GetLatticeManager()->GetLattice(aStep.GetPostStepPoint()->GetPhysicalVolume())
	   << G4endl;
  }
  
  //First, check the different materials involved. If one of them does not have
  //a lattice, then something upstream has gone wrong. Scream. Note that here,
  //we need to use the G4LatticeManager rather than the TrackUtils version of
  //GetLattice() because the trackUtils version somehow chooses a point that's
  //inside the current (i.e. not the post-step) volume.
  G4LatticePhysical * latNear = G4LatticeManager::GetLatticeManager()->GetLattice(aStep.GetPreStepPoint()->GetPhysicalVolume());
  G4LatticePhysical * latFar = G4LatticeManager::GetLatticeManager()->GetLattice(aStep.GetPostStepPoint()->GetPhysicalVolume());
  if (!latNear || !latFar) {
    G4ExceptionDescription msg;
    msg << "Expecting to do phonon transmission at interface but one or more "
	<< "lattices splitting the boundary cannot be found.";
    G4Exception("G4CMPPhononBoundaryProcess::DoTransmission",
		"PhononBoundary001",FatalException,msg);
  }

  //SIMPLEST POSSIBLE IMPLEMENTATION -- PHYSICS IS NOT NECESSARILY RIGHT
  //Now get track info at this point: the k vector, mode, etc.
  auto trackInfo = G4CMP::GetTrackInfo<G4CMPPhononTrackInfo>(aTrack);
  G4ThreeVector waveVector = trackInfo->k();
  G4int mode = GetPolarization(aStep.GetTrack());
  
  // Check whether step has proper boundary-stopped geometry
  G4ThreeVector surfacePoint;
  if (!CheckStepBoundary(aStep, surfacePoint)) {
    if (verboseLevel > 5) {
      G4cout << "DT Function Point B | checking step boundary failed in "
	     << "DoTransmission" << G4endl;
    }
    if (verboseLevel>2)
      G4cout << " Boundary point moved to " << surfacePoint << G4endl;

    aParticleChange.ProposePosition(surfacePoint);	// IS THIS CORRECT?!?
  }

  
  //Since the lattice hasn't changed yet, change it here. (This also happens
  //at the MFP calc point at the beginning of the next step, but it's nice to
  //have it here so we can use the new lattice info to help figure out vdir,
  //etc.). Also need to update the lattice for anharmonic decay, since that's
  //a data member here and isn't getting updated elsewhere.)
  this->SetLattice(G4LatticeManager::GetLatticeManager()->GetLattice(aStep.GetPostStepPoint()->GetPhysicalVolume()));
  anharmonicDecay->SetLattice(G4LatticeManager::GetLatticeManager()->GetLattice(aStep.GetPostStepPoint()->GetPhysicalVolume()));
  UpdateSCAfterLatticeChange();

  //Debugging
  if( verboseLevel > 5 ){    
    G4cout << "DT Function Point C | the lattice at the end of doTransmission: "
	   << theLattice << G4endl;
  }


  //Need to rotate to local coordinates. Since the current touchable at this
  //point is still the pre-transmission volume, we will need to use the
  //post-step touchable for rotations
  const G4VTouchable * nextVolTouchable
    = aStep.GetPostStepPoint()->GetTouchable();

  //Inbound and outbound wavevectors are the same (waveVector) by aphysical
  //fiat at the moment -- Caitlyn will change this
  trackInfo->SetWaveVector(waveVector); 
  G4CMP::RotateToLocalDirection(nextVolTouchable,waveVector);  
  G4ThreeVector vdir = theLattice->MapKtoVDir(mode, waveVector);
  G4double v = theLattice->MapKtoV(mode, waveVector);
  G4CMP::RotateToGlobalDirection(nextVolTouchable,vdir);  
  aParticleChange.ProposeVelocity(v);
  aParticleChange.ProposeMomentumDirection(vdir);

  //Rotate the wavevector back to global dimensions so we can have it back for
  //further calculations
  G4CMP::RotateToGlobalDirection(aStep.GetPostStepPoint()->GetTouchable(),waveVector);  

  //Get the surface normal, compute a generalizedSurfaceNormal from it (to
  //indicate the direction pointing with the incident phonon),
  //and then check the direction of the outgoing phonon with respect to it
  G4ThreeVector incMomDir = aStep.GetPreStepPoint()->GetMomentumDirection();
  G4ThreeVector generalizedSurfNorm =
    G4CMP::GetGeneralizedSurfaceNormal(G4CMP::GetSurfaceNormal(aStep),
				       incMomDir);
  
  //If the direction is not "outward" with respect to the direction from which
  //the phonon came, then we have an issue.
  if (G4CMP::PhononVelocityIsOutward(theLattice,mode,waveVector,generalizedSurfNorm,nextVolTouchable,surfacePoint) == false) {

    //For now, let's do something unphysical and kludgey: pick new k-vectors
    //within a small angle of the original one to be transmitted, recompute the
    //direction, and retry the check. This will have to be fixed by some
    //combination of Caitlyn and others
    int nTries = 1000;
    int currentTry = 0;
    G4ThreeVector newAttemptWaveVector;
    G4bool successfulTransmission = false;
    for( int iT = 0; iT < nTries; ++iT ){

      //Make the length approximately 40% of the wavevector.
      //This seems egregious...
      G4ThreeVector adjustmentVector =
	G4RandomDirection() * waveVector.mag() * 0.4; 

      //Do some math to at least keep the length of the vector the same
      G4double oldWaveVectorMag = waveVector.mag();
      G4double newAttemptWaveVectorMag = (waveVector+adjustmentVector).mag();     
      newAttemptWaveVector = (waveVector + adjustmentVector) *
	oldWaveVectorMag / newAttemptWaveVectorMag;

      //If the phonon is now pointing outward, then we can store particle
      //change info and overwrite vDir and v
      if (G4CMP::PhononVelocityIsOutward(theLattice,mode,newAttemptWaveVector,
					 generalizedSurfNorm,nextVolTouchable,
					 surfacePoint) == true) {

	//Sanity check
	if (verboseLevel > 5) {    
	  G4cout << "DT Function Point D | in success trigger of "
		 << "phononVelIsOutward loop, newAttemptWaveVector: "
		 << newAttemptWaveVector << G4endl;
	}
	
	//Put the newAttemptWaveVector into a local direction within the new
	//touchable's geometry, then do mappings
	G4CMP::RotateToLocalDirection(nextVolTouchable,newAttemptWaveVector); 	
	vdir = theLattice->MapKtoVDir(mode, newAttemptWaveVector);
	v = theLattice->MapKtoV(mode, newAttemptWaveVector);

	if (verboseLevel > 5) {    
	  G4cout << "DT Function Point E | in success trigger of "
		 << "phononVelIsOutward loop, newAttemptWaveVector rotated "
		 << "into local: " << newAttemptWaveVector << G4endl;
	}

	//Put the newAttemptWaveVector and the vdir back into a global direction
	G4CMP::RotateToGlobalDirection(nextVolTouchable,newAttemptWaveVector);
	G4CMP::RotateToGlobalDirection(nextVolTouchable,vdir);

	if (verboseLevel > 5) {    
	  G4cout << "DT Function Point F | in success trigger of "
		 << "phononVelIsOutward loop, newAttemptWaveVector rotated "
		 << "back into global: " << newAttemptWaveVector << G4endl;
	}
	
	//Now fill out quantities
	//Inbound and outbound wavevectors are the same (waveVector) by
	//aphysical fiat at the moment.
	trackInfo->SetWaveVector(newAttemptWaveVector); 
	aParticleChange.ProposeVelocity(v);
	aParticleChange.ProposeMomentumDirection(vdir);
	successfulTransmission = true;
	break;
      }
    }

    if (verboseLevel > 5) {
      G4cout << "DT Function Point G | Momentum direction proposed at end of "
	     << "doTransmission: " << vdir << G4endl;
    }
    
    //If after our loop we still run into issues, we have to break and kill
    //the track
    if (!successfulTransmission) {
      G4String msg = G4PhononPolarization::Name(mode) +
	" transmission failed even after repeated tries to nudge direction.";

      if (verboseLevel > 5) {
	G4cout << "DT Function Point H | Phonon incident waveVector: "
	       << waveVector.unit() << ", generalized surfNorm: "
	       << generalizedSurfNorm << ", attempted transmitted k-vector "
	       << waveVector.unit() << ", and attempted transmitted vDir: "
	       << vdir << G4endl;
      }
      G4Exception((GetProcessName()+"::DoTransmission").c_str(), "Boundary011",
		  JustWarning, msg.c_str());
      DoSimpleKill(aTrack, aStep, aParticleChange);
    }
    return;
  }    
}

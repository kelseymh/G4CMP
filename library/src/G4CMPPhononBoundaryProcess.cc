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

#include "G4CMPPhononBoundaryProcess.hh"
#include "G4CMPAnharmonicDecay.hh"
#include "G4CMPConfigManager.hh"
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

G4double G4CMPPhononBoundaryProcess::GetMeanFreePath(const G4Track& /*aTrack*/,
                                             G4double /*prevStepLength*/,
                                             G4ForceCondition* condition) {
  *condition = Forced;
  return DBL_MAX;
}


// Process action

G4VParticleChange*
G4CMPPhononBoundaryProcess::PostStepDoIt(const G4Track& aTrack,
                                         const G4Step& aStep) {
  // NOTE:  G4VProcess::SetVerboseLevel is not virtual!  Can't overlaod it
  G4CMPBoundaryUtils::SetVerboseLevel(verboseLevel);

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

  ApplyBoundaryAction(aTrack, aStep, phParticleChange);
  
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

    G4ThreeVector vec1 = G4CMP::GetLambertianVector(theLattice, surfNorm, mode,
                                                    surfacePoint);
    G4ThreeVector vec2 = G4CMP::GetLambertianVector(theLattice, surfNorm, mode,
                                                    surfacePoint);

    sec1->SetMomentumDirection(vec1);
    sec2->SetMomentumDirection(vec2);

    return;
  } else if (random < downconversionProb + specProb) {
    reflectedKDir = GetSpecularVector(waveVector, surfNorm, mode, surfacePoint); // Modify surfacePoint & surfNorm in place
    refltype = "specular";
  } else {
    reflectedKDir = G4CMP::GetLambertianVector(theLattice, surfNorm, mode,
                                               surfacePoint);
    refltype = "diffuse";
  }

  // Update trackInfo wavevector and particleChange's group velocity and momentum direction
  // reflectedKDir is in global coordinates here - no conversion needed
  FillParticleChange(particleChange, aTrack, reflectedKDir);
  const G4ThreeVector vdir = *particleChange.GetMomentumDirection();

  // If displacement occured: update the particle change and navigator for volume assignment
  if (refltype == "specular" && *particleChange.GetPosition() != surfacePoint) {
    FillParticleChange(phParticleChange, aStep, surfacePoint);
    UpdateNavigatorVolume(aStep, surfacePoint, vdir);
  }

  if (verboseLevel>2) {
    G4cout << "New surface position " << *particleChange.GetPosition()/mm 
	   << " mm\n New wavevector direction " << reflectedKDir
	   << "\n New momentum direction " << vdir << G4endl;
  }

  // If reflection failed, report problem and kill the track
  if (!G4CMP::PhononVelocityIsInward(theLattice,mode,reflectedKDir,surfNorm, surfacePoint)) {
    G4String msg = G4PhononPolarization::Name(mode) + " " + refltype + " reflection failed";

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
                  G4ThreeVector& surfacePoint) {
  // Specular reflecton should reverse momentum along normal
  G4ThreeVector reflectedKDir = waveVector.unit();
  G4double kPerp = reflectedKDir * surfNorm;		// Dot product between k and norm
  (reflectedKDir -= 2.*kPerp*surfNorm).setMag(1.);	// Reflect against normal

  if (G4CMP::PhononVelocityIsInward(theLattice,mode,reflectedKDir,surfNorm,
                                    surfacePoint))
    return reflectedKDir;

  // Reflection didn't work as expected, need to correct:
  // If the reflected wave vector cannot propagate in the bulk
  // (i.e., the reflected k⃗ has an associated v⃗g which is not inwardly directed.)
  // That surface wave will propagate until it reaches a point
  // where the wave vector has an inwardly directed v⃗g.
  RotateToLocalDirection(reflectedKDir);
  G4ThreeVector newNorm = surfNorm;
  RotateToLocalDirection(newNorm);
  
  // Initialize solid object and utilities
  G4VSolid* solid = GetCurrentVolume()->GetLogicalVolume()->GetSolid();
  G4CMPSolidUtils solidUtils(solid, verboseLevel, GetProcessName());

  G4ThreeVector stepLocalPos = GetLocalPosition(surfacePoint);
  G4ThreeVector oldNorm = newNorm;
  G4ThreeVector oldstepLocalPos = stepLocalPos;

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
      // Reflect kTan against the edge - rotates & modifies kTan; modifies newNorm
      solidUtils.ReflectAgainstEdge(kTan, stepLocalPos, newNorm);
    } else {
      // Rotate kTan to new position
      axis = newNorm.cross(kTan).unit();
      phi = oldNorm.azimAngle(newNorm, axis);
      kTan = kTan.rotate(axis, phi);
    }

    // Get perpendicular component of reflected k w/ new norm
    // (negative implied in kPerpMag for inward pointing)
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

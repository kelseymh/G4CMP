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

#include "G4CMPPhononBoundaryProcess.hh"
#include "G4CMPAnharmonicDecay.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPGeometryUtils.hh"
#include "G4CMPPhononTrackInfo.hh"
#include "G4CMPSurfaceProperty.hh"
#include "G4CMPTrackUtils.hh"
#include "G4CMPUtils.hh"
#include "G4CMPConfigManager.hh"
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

  aParticleChange.Initialize(aTrack);

  // Handle Boundary -> Boundary steps that "ignore" the bulk of the detector
  // e.g. DetectorPV -> DetectorPV instead of DetectorPV -> ZipPV -> DetectorPV
  if (BoundaryToBoundaryStep(aStep)) {
    GetBoundingVolumes(aStep);
    SetPrePV(GetCurrentVolume());
    GetSurfaceProperty(aStep);
  } else if (!IsGoodBoundary(aStep)) {
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }

  if (verboseLevel>1) {
    G4cout << GetProcessName() << "::PostStepDoIt "
           << "Event " << G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID()
           << " Track " << aTrack.GetTrackID() << " Step " << aTrack.GetCurrentStepNumber()
           << G4endl;
  }

  if (verboseLevel>2) {
    G4cout << " K direction: " << GetLocalWaveVector(aTrack).unit()
           << "\n P direction: " << aTrack.GetMomentumDirection() << G4endl;
  }

  ApplyBoundaryAction(aTrack, aStep, aParticleChange);

  ClearNumberOfInteractionLengthLeft();		// All processes should do this!
  G4cout << "---------------------------------------------------------------" << G4endl;
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

G4bool G4CMPPhononBoundaryProcess::BoundaryToBoundaryStep(const G4Step& aStep) {
  /*G4StepPoint* preStep = aStep.GetPreStepPoint();
  G4StepPoint* postStep = aStep.GetPostStepPoint();
  return (preStep->GetStepStatus() == fGeomBoundary &&
          postStep->GetStepStatus() == fGeomBoundary &&
          preStep->GetPhysicalVolume() != GetCurrentVolume() &&
          preStep->GetPhysicalVolume() == postStep->GetPhysicalVolume());*/
  return false;
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
  G4ThreeVector surfacePoint = aStep.GetPostStepPoint()->GetPosition();
  G4ThreeVector surfNorm = G4CMP::GetSurfaceNormal(aStep);

  if (verboseLevel>2) {
    G4cout << "\n Old wavevector direction " << waveVector.unit()
	   << "\n Old momentum direction   " << aTrack.GetMomentumDirection()
	   << G4endl;
  }

  // Check whether step has proper boundary-stopped geometry
  if (!BoundaryToBoundaryStep(aStep) && !CheckStepBoundary(aStep, surfacePoint)) {
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
    reflectedKDir = GetReflectedVector(waveVector, surfNorm, mode, surfacePoint); // Modify surfacePoint & surfNorm in place
    refltype = "specular";
  } else {
    reflectedKDir = GetLambertianVector(surfNorm, mode);
    refltype = "diffuse";
  }

  // Update trackInfo wavevector and particleChange's group velocity and momentum direction
  // reflectedKDir is in global coordinates here - no conversion needed
  FillParticleChange(particleChange, aTrack, reflectedKDir);
  const G4ThreeVector vdir = *particleChange.GetMomentumDirection();

  // If displacement occured: update the surface position and navigator for volume assignment
  if (*particleChange.GetPosition() != surfacePoint) {
    G4cout << "CHANGING POSITION AND UPDATING NAV" << G4endl;
    particleChange.ProposePosition(surfacePoint);
    G4Navigator* navigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
    navigator->LocateGlobalPointAndSetup(surfacePoint, &vdir, true, false);
  }

  if (verboseLevel>2) {
    G4cout << "\n New surface position " << *particleChange.GetPosition()
     << "\n New wavevector direction " << reflectedKDir
	   << "\n New momentum direction " << vdir << G4endl;
  }

  // If reflection failed, report problem and kill the track
  if (!G4CMP::PhononVelocityIsInward(theLattice,mode,reflectedKDir,surfNorm)) {
    G4Exception((GetProcessName()+"::DoReflection").c_str(), "Boundary010",
		JustWarning, ("Phonon "+refltype+" reflection failed"+"\nPhonon mode at time of death: "+G4PhononPolarization::Label(mode)).c_str());
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
    G4VSolid* solid = GetCurrentVolume()->GetLogicalVolume()->GetSolid();

    EInside place = solid->Inside(stepLocal);
    G4cout << " After trial step, " << (place==kInside ? "inside"
					: place==kOutside ? "OUTSIDE"
					: "on surface") << G4endl;
  }
}


// Generate specular reflection corrected for momentum dispersion

G4ThreeVector G4CMPPhononBoundaryProcess::
GetReflectedVector(const G4ThreeVector& waveVector,
		               G4ThreeVector& surfNorm, G4int mode,
		               G4ThreeVector& surfacePoint) {
  // Specular reflecton should reverse momentum along normal
  G4ThreeVector reflectedKDir = waveVector.unit();
  G4double kPerp = reflectedKDir * surfNorm;		// Dot product between k and norm
  (reflectedKDir -= 2.*kPerp*surfNorm).setMag(1.);	// Reflect against normal

  if (G4CMP::PhononVelocityIsInward(theLattice,mode,reflectedKDir,surfNorm))
    return reflectedKDir;

  // Reflection didn't work as expected, need to correct:
  // If the reflected wave vector cannot propagate in the bulk
  // (i.e., the reflected k⃗ has an associated v⃗g which is not inwardly directed.)
  // That surface wave will propagate until it reaches a point
  // where the wave vector has an inwardly directed v⃗g.
  RotateToLocalDirection(reflectedKDir);
  G4ThreeVector newNorm = surfNorm;
  RotateToLocalDirection(newNorm);
  
  G4ThreeVector stepLocalPos = GetLocalPosition(surfacePoint);
  G4VSolid* solid = GetCurrentVolume()->GetLogicalVolume()->GetSolid();
  G4ThreeVector oldNorm = newNorm;
  // Find the distance from point to surface along norm (- means inward)
  G4double surfAdjust = solid->DistanceToIn(stepLocalPos, -newNorm);
  G4double kPerpMag = reflectedKDir.dot(newNorm);

  G4ThreeVector kPerpV = kPerpMag * newNorm;		// Negative implied in kPerpMag for inward pointing
  G4ThreeVector kTan = reflectedKDir - kPerpV;		// Get kTan: reflectedKDir = kPerpV + kTan
  G4double kTanMag = kTan.mag();
  EInside isIn = solid->Inside(stepLocalPos);

  G4int nAttempts = 0;

  // debugging only DELETE
  G4ThreeVector oldkTan = kTan;
  G4ThreeVector oldkPerpV = kPerpV;
  G4ThreeVector oldstepLocalPos = stepLocalPos;

  // Initialize stepSize for _this_ solid object
  G4CMPConfigManager* config = G4CMPConfigManager::Instance();
  stepSize = config->GetPhononSurfStepSize();
  // Set default stepSize based on solid bounding limits
  if (stepSize == 0) {
    G4ThreeVector pmin(0,0,0);
    G4ThreeVector pmax(0,0,0);
    solid->BoundingLimits(pmin, pmax);
    stepSize = (pmax - pmin).mag() / nStepLimit;
  }

  // FIXME: Need defined units
  if (verboseLevel>3) {
    G4cout << "GetReflectedVector:beforeLoop -> "
      << ", stepLocalPos = " << stepLocalPos
      << ", kPerpMag (newNorm dot reflectedKDir) = " << kPerpMag
      << ", newNorm = " << newNorm
      << ", reflectedKDir = " << reflectedKDir
      << ", kPerpV (kPerpMag * newNorm) = " << kPerpV
      << ", kTan (reflectedKDir - kPerpV) = " << kTan
      << ", surfaceStepSize = " << G4BestUnit(stepSize, "Length")
      << ", nStepLimit = " << nStepLimit << G4endl;
  }

  // Assumes everything is in Global. Just add the GetGlobal in the loop conditions.
  while (!G4CMP::PhononVelocityIsInward(theLattice, mode,
   GetGlobalDirection(reflectedKDir), GetGlobalDirection(newNorm))
	 && nAttempts++ < nStepLimit) {
    // Save previous loop values
    oldstepLocalPos = stepLocalPos;
    oldNorm = newNorm;
    // debugging only - DELETE
    oldkTan = kTan;
    oldkPerpV = kPerpV;

    // Step along kTan direction - this point is now outside the detector
    stepLocalPos += stepSize * kTan.unit();

    // Get the local normal at the new surface point
    newNorm = solid->SurfaceNormal(stepLocalPos);

    // Adjust stepLocalPos back to surface of detector
    surfAdjust = solid->DistanceToIn(stepLocalPos, -newNorm);
    stepLocalPos -= surfAdjust * newNorm;
    isIn = solid->Inside(stepLocalPos);

    // Large normal changes and not being on surface after initial adjustment
    // indicates we are approaching an edge
    if (isIn != kSurface || newNorm * oldNorm <= 0) {
      // Reset stepLocalPos and newNorm to last valid surface point
      stepLocalPos = oldstepLocalPos;
      newNorm = oldNorm;
      // Modify stepLocalPos in place to edge position
      AdjustToEdgePosition(solid, kTan, stepLocalPos);
      // Reflect kTan against the edge - rotates & modifies kTan; modifies newNorm
      ReflectAgainstEdge(solid, kTan, stepLocalPos, newNorm);
    }
    else {
      // "Rotate" kTan to new position
      (kTan -= newNorm * (kTan * newNorm)).setMag(kTanMag);
    }

    // Get perpendicular component of reflected k w/ new norm (negative implied in kPerpMag for inward pointing)
    kPerpV = kPerpMag * newNorm;

    // Calculate new reflectedKDir (kTan + kPerpV)
    reflectedKDir = kTan + kPerpV;

    // Debugging: Can be removed?
    G4ThreeVector vDir = theLattice->MapKtoVDir(mode, reflectedKDir);

    // FIXME: Need defined units
    if (verboseLevel>3) {
      G4cout << " "
       << "GetReflectedVector:insideLoop -> "
       << "attempts = " << nAttempts
       << ", oldstepLocalPos = " << oldstepLocalPos
       << ", surfAdjust = " << surfAdjust
       << ", stepLocalPos = " << stepLocalPos
       << ", oldkPerpV = " << oldkPerpV
       << ", oldkTan = " << oldkTan
       << ", kPerpV (kPerpMag * newNorm) = " << kPerpV
       << ", kPerpMag = " << kPerpMag
       << ", newNorm = " << newNorm
       << ", oldNorm = " << oldNorm
       << ", kTan = " << kTan
       << ", reflectedKDir (kTan + kPerpV) = " << reflectedKDir
       << ", Phonon mode = " << G4PhononPolarization::Label(mode)
       << ", New group velocity: " << vDir << G4endl;
    }
  }

  // Restore global coordinates to reflectedKDir & newNorm
  RotateToGlobalDirection(reflectedKDir);
  RotateToGlobalDirection(newNorm);

  if (!G4CMP::PhononVelocityIsInward(theLattice, mode, reflectedKDir, newNorm)) {
    G4cout << (GetProcessName()+"::GetReflectedVector").c_str()
      << ": Phonon displacement failed after " << nAttempts - 1 
      << " attempts. Doing diffuse reflection at surface point: " << surfacePoint << G4endl;

    // reflectedKDir and stepLocalPos will be in global coordinates
    reflectedKDir = GetLambertianVector(surfNorm, mode);
    stepLocalPos = surfacePoint;
  }
  else{ 
    // Restore global coordinates to stepLocalPos
    RotateToGlobalPosition(stepLocalPos);
    // Update surfNorm to new point's normal
    surfNorm = newNorm;
  }

  if (verboseLevel>2) {
    G4cout << (GetProcessName()+"::GetReflectedVector").c_str()
      << ": nAttempts = " << nAttempts
      << ", waveVector = " << waveVector
      << ", reflectedKDir = " << reflectedKDir
      << ", initialGlobalPostion = " << surfacePoint
      << ", finalGlobalPosition = " << stepLocalPos << G4endl;
  }

  surfacePoint = stepLocalPos;
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


// Find the direction for minimum distance to surface

void G4CMPPhononBoundaryProcess::
OptimizeSurfaceAdjustAngle(const G4VSolid* solid,
                           const G4ThreeVector& stepLocalPos, G4double& theta0,
                           G4double& phi0, const G4int angOption,
                           const G4double minDist) const {
  // Constants used for Golden Section searches
  G4double const tolerance = 1e-12*rad;
  G4double const gRatio = (std::sqrt(5) + 1) / 2;
  EInside const isIn = solid->Inside(stepLocalPos);

  // Define search quantities for finding optimal angle
  G4double a = 0;
  G4double b = (angOption == 0) ? pi : twopi;
  G4double x1 = a + (b - a) / gRatio;
  G4double x2 = b - (b - a) / gRatio;
  G4ThreeVector dir1(0,0,0);
  G4ThreeVector dir2(0,0,0);
  G4double coeffX = 0, coeffY = 0, coeffZ = 0;

  // Initial search directions
  if (angOption == 0) {
    // Optimize theta
    coeffX = minDist*cos(phi0);
    coeffY = minDist*sin(phi0);
    coeffZ = minDist;
    dir1.set(coeffX*sin(x1), coeffY*sin(x1), coeffZ*cos(x1));
    dir2.set(coeffX*sin(x2), coeffY*sin(x2), coeffZ*cos(x2));
  }
  else {
    // Optimize phi
    coeffX = minDist*sin(theta0);
    coeffY = minDist*sin(theta0);
    coeffZ = minDist*cos(theta0);
    dir1.set(coeffX*cos(x1), coeffY*sin(x1), coeffZ);
    dir2.set(coeffX*cos(x2), coeffY*sin(x2), coeffZ);
  }

  // Distances for each direction
  G4double dist1 = (isIn == kInside) ? solid->DistanceToOut(stepLocalPos + dir1)
                                     : solid->DistanceToIn(stepLocalPos + dir1);
  G4double dist2 = (isIn == kInside) ? solid->DistanceToOut(stepLocalPos + dir2)
                                     : solid->DistanceToIn(stepLocalPos + dir2);

  // Golden section search for optimizing angle
  while (b - a > tolerance) {
    // Adjust bounds
    if (dist1 < dist2) {
      // Shift up
      a = x2;
      x2 = x1;
      dist2 = dist1;
      x1 = a + (b - a) / gRatio;
      if (angOption == 0) {
        // Optimize theta
        dir1.set(coeffX*sin(x1), coeffY*sin(x1), coeffZ*cos(x1));
      }
      else {
        // Optimize phi
        dir1.set(coeffX*cos(x1), coeffY*sin(x1), coeffZ);
      }
      dist1 = (isIn == kInside) ? solid->DistanceToOut(stepLocalPos + dir1)
                                : solid->DistanceToIn(stepLocalPos + dir1);
    }
    else {
      // Shift down
      b = x1;
      x1 = x2;
      dist1 = dist2;
      x2 = b - (b - a) / gRatio;
      if (angOption == 0) {
        // Optimize theta
        dir2.set(coeffX*sin(x2), coeffY*sin(x2), coeffZ*cos(x2));
      }
      else {
        // Optimize phi
        dir2.set(coeffX*cos(x2), coeffY*sin(x2), coeffZ);
      }
      dist2 = (isIn == kInside) ? solid->DistanceToOut(stepLocalPos + dir2)
                                : solid->DistanceToIn(stepLocalPos + dir2);
    }
  }
  G4double optimalAng = (x1 + x2) / 2;

  // Adjust angle in place
  if (angOption == 0) { theta0 = optimalAng; }
  else { phi0 = optimalAng; }
}


// Adjust to surface position closest to stepLocalPos without SurfaceNormal

void G4CMPPhononBoundaryProcess::
AdjustToClosestSurfacePoint(const G4VSolid* solid,
                            G4ThreeVector& stepLocalPos) const {
  // Determine where you are in respect to solid
  EInside isIn = solid->Inside(stepLocalPos);
  if (isIn == kSurface) return;

  // Angles to be adjusted in place by OptimizeSurfaceAdjustAngle
  // Best initial conditions will be the surface normal
  G4ThreeVector surfNorm = solid->SurfaceNormal(stepLocalPos);
  G4double bestTheta = surfNorm.theta();
  G4double bestPhi = surfNorm.phi();

  G4double minDist = (isIn == kInside) ? solid->DistanceToOut(stepLocalPos)
                                       : solid->DistanceToIn(stepLocalPos);

  // Try to optimize phi first
  OptimizeSurfaceAdjustAngle(solid, stepLocalPos, bestTheta, bestPhi, 1, minDist);
  OptimizeSurfaceAdjustAngle(solid, stepLocalPos, bestTheta, bestPhi, 0, minDist);

  G4ThreeVector optDir(minDist*sin(bestTheta)*cos(bestPhi),
                       minDist*sin(bestTheta)*sin(bestPhi),
                       minDist*cos(bestTheta));

  // Set stepLocalPos and exit if a surface point was found
  if (solid->Inside(stepLocalPos + optDir) == kSurface) {
    stepLocalPos += optDir;
    return;
  }

  // Retry but optimize theta first
  bestTheta = surfNorm.theta();
  bestPhi = surfNorm.phi();

  OptimizeSurfaceAdjustAngle(solid, stepLocalPos, bestTheta, bestPhi, 0, minDist);
  OptimizeSurfaceAdjustAngle(solid, stepLocalPos, bestTheta, bestPhi, 1, minDist);

  optDir.set(minDist*sin(bestTheta)*cos(bestPhi),
             minDist*sin(bestTheta)*sin(bestPhi),
             minDist*cos(bestTheta));

  // Only return valid positions on surface
  if (solid->Inside(stepLocalPos + optDir) == kSurface) {
    stepLocalPos += optDir;
  }
  else {
    stepLocalPos.set(kInfinity,kInfinity,kInfinity);
  }
}


// Do a binary search to find the closest point toward the edge along kTan
// Modifies stepLocalPos in place

void G4CMPPhononBoundaryProcess::
AdjustToEdgePosition(const G4VSolid* solid, const G4ThreeVector& kTan,
                     G4ThreeVector& stepLocalPos) const {
  EInside isIn = solid->Inside(stepLocalPos);
  G4double low = 0.0*um;
  G4double high = stepSize;
  G4double mid = 0;
  G4ThreeVector originalPos = stepLocalPos;
  G4double tolerance = solid->GetTolerance();

  // Binary search to bring surface point to edge
  while (high - low > tolerance || isIn != kSurface) {
    mid = 0.5 * (low + high);
    stepLocalPos = originalPos + mid * kTan;
    // Modify stepLocalPos in place
    AdjustToClosestSurfacePoint(solid, stepLocalPos);

    isIn = solid->Inside(stepLocalPos);

    if (isIn == kSurface) {
      low = mid; // Move out
    }
    else {
      high = mid; // Move in
    }
  }

  if (verboseLevel>3) {
    G4double remDist = (solid->DistanceToIn(stepLocalPos) >=
                        solid->DistanceToOut(stepLocalPos))
                      ? solid->DistanceToIn(stepLocalPos)
                      : solid->DistanceToOut(stepLocalPos);
    G4cout << (GetProcessName()+"::AdjustToEdgePosition").c_str()
      << ": initialPos = " << originalPos
      << ", kTan = " << kTan
      << ", finalPos = " << stepLocalPos
      << ", remainingDist = " << remDist << G4endl;
  }
}


// Reflect kTan against an edge
// Modifies kTan and newNorm in place

void G4CMPPhononBoundaryProcess::
ReflectAgainstEdge(const G4VSolid* solid, G4ThreeVector& kTan,
                   const G4ThreeVector& stepLocalPos, G4ThreeVector& newNorm) const {
  // Get normal of both surfaces at this edge point by stepping with normAdjust
  G4double normAdjust = 1*nm;
  G4double kTanMag = kTan.mag();
  G4ThreeVector norm1 = solid->SurfaceNormal(stepLocalPos);
  G4ThreeVector norm2 = solid->SurfaceNormal(stepLocalPos - normAdjust*norm1);
  // Try to fix norm1 if we didn't get the normals for a proper reflection
  norm1 = (norm1 * norm2 > 0) ? solid->SurfaceNormal(stepLocalPos - normAdjust*norm2) : norm1;

  G4ThreeVector edgeVec(0,0,0);
  G4ThreeVector refNorm(0,0,0);

  // Only do reflection if at an edge
  if (norm1 * norm2 <= 0) {
    newNorm = (newNorm*norm1 > newNorm*norm2) ? norm1 : norm2;
    // Rotate kTan to be orthogonal to one normal
    (kTan -= newNorm * (kTan * newNorm)).setMag(kTanMag);

    // Get the edge vector
    edgeVec = norm1.cross(norm2).unit();

    // Find the normal for the reflection "surface"
    refNorm = (edgeVec.cross(newNorm)).unit();
    if (refNorm * kTan < 0) refNorm *= -1;

    // Reflect kTan against reflection "surface"
    kTan -= 2*((kTan * refNorm) * refNorm);
  }

  if (verboseLevel>3) {
    G4cout << (GetProcessName()+"::ReflectAgainstEdge").c_str()
      << ": stepLocalPos = " << stepLocalPos
      << ", kTan_0 = " << kTan - 2*((kTan * -refNorm) * -refNorm)
      << ", edgeVector = " << edgeVec
      << ", refNorm = " << refNorm
      << ", norm1 = " << norm1
      << ", norm2 = " << norm2
      << ", kTan_f = " << kTan << G4endl;
  }
}

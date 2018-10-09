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

#include "G4CMPPhononBoundaryProcess.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPGeometryUtils.hh"
#include "G4CMPPhononTrackInfo.hh"
#include "G4CMPSurfaceProperty.hh"
#include "G4CMPTrackUtils.hh"
#include "G4CMPUtils.hh"
#include "G4ExceptionSeverity.hh"
#include "G4GeometryTolerance.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4Navigator.hh"
#include "G4ParallelWorldProcess.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4TransportationManager.hh"
#include "G4VParticleChange.hh"
#include "G4LogicalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4CMPAnharmonicDecay.hh"
#include "Randomize.hh"


G4CMPPhononBoundaryProcess::G4CMPPhononBoundaryProcess(const G4String& aName)
  : G4VPhononProcess(aName, fPhononReflection), G4CMPBoundaryUtils(this) {;}

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
	   << "\n k_perp " << k*G4CMP::GetSurfaceNormal(aStep)
	   << " vs. absMinK " << absMinK << G4endl;
  }

  return (G4CMPBoundaryUtils::AbsorbTrack(aTrack,aStep) &&
    k*G4CMP::GetSurfaceNormal(aStep) > absMinK);
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
    G4cout << " Old wavevector direction " << waveVector.unit()
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

  G4double Eoverh = GetKineticEnergy(aTrack)/h_Planck;
  G4double EoverhGHz = Eoverh * ns;

  G4double specProb = BoundarySpecularProb(EoverhGHz);
  G4double diffuseProb = BoundaryLambertianProb(EoverhGHz);
  G4double downconversionProb = BoundaryAnharmonicProb(EoverGHz);

  // Empirical functions may lead to non normalised probabilities.
  // Normalise here.

  G4double norm = specProb + diffuseProb + downconversionProb;

  specProb /= norm;
  diffuseProb /= norm;
  downconversionProb /= norm;

  G4ThreeVector reflectedKDir;

  G4double random = G4G4UniformRand();

  if (random < downconversionProb) {
    /* Do Downconversion */
    G4G4CMPAnharmonicDecay::DoDecay(aTrack, aStep, particleChange);
    G4Track* sec1 = particleChange.GetSecondary(0);
    G4Track* sec2 = particleChange.GetSecondary(1);

    G4ThreeVector vec1 = GetLambertianVector(surfNorm);
    G4ThreeVector vec2 = GetLambertianVector(surfNorm);

    sec1->SetMomentumDirection(vec1);
    sec2->SetMomentumDirection(vec2);
    return;

  } else if (random < downconversionProb + specProb) {
    // Specular reflecton reverses momentum along normal
    reflectedKDir = waveVector.unit();
    G4double kPerp = reflectedKDir * surfNorm;
    reflectedKDir -= 2.*kPerp * surfNorm;
  } else {
    reflectedKDir = GetLambertianVector(surfNorm);
  }

  // If reflection failed, report problem and kill the track
  if (!G4CMP::PhononVelocityIsInward(theLattice,mode,reflectedKDir,surfNorm)) {
    G4Exception((GetProcessName()+"::DoReflection").c_str(), "Boundary010",
		JustWarning, "Phonon reflection failed");
    DoSimpleKill(aTrack, aStep, aParticleChange);
    return;
  }

  G4ThreeVector vdir = theLattice->MapKtoVDir(mode, reflectedKDir);
  G4double v = theLattice->MapKtoV(mode, reflectedKDir);

  if (verboseLevel>2) {
    G4cout << " New wavevector direction " << reflectedKDir
	   << "\n New momentum direction   " << vdir << G4endl;
  }

  // SANITY CHECK:  Project a 1 um step in the new direction, see if it
  // is still in the correct (pre-step) volume.

  if (verboseLevel>2) {
    G4ThreeVector stepPos = surfacePoint + .1*mm * vdir;

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

G4double G4CMPPhononBoundaryProcess::BoundaryAnharmonicProb(const G4double f_GHz) {

  // 520 GHz is where probabilities are undefined, just use previous
  // probabilities assuming no downconversion.
  if (f_GHz > 520) return 0.0;
  return 1.51e-14 * (f_GHz * f_GHz * f_GHz * f_GHz * f_GHz);
}

G4double G4CMPPhononBoundaryProcess::BoundarySpecularProb(const G4double f_GHz) {
  // 350 GHz unphysical cutoff, probably should define this as some variable
  if (f_GHz > 520) return GetMaterialProperty("specProb");
  if (f_GHz >= 350) return 1 - BoundaryLambertianProb(350) - BoundaryAnharmonicProb(f_GHz);
  return 2.9e-13 * (f_GHz * f_GHz * f_GHz * f_GHz) +
         3.1e-9 * (f_GHz * f_GHz * f_GHz) -
         3.21e-6 * (f_GHz * f_GHz) -
         2.03e-4 * f_GHz +
         0.928;
}

G4double G4CMPPhononBoundaryProcess::BoundaryLambertianProb(const G4double f_GHz) {
  if (f_GHz > 520) return 1 - GetMaterialProperty("specProb");
  if (f_GHz >= 350) f_GHz = 350;
  return -2.98e-11 * (f_GHz * f_GHz * f_GHz * f_GHz) +
         1.71e-8 * (f_GHz * f_GHz * f_GHz) -
         2.47e-6 * (f_GHz * f_GHz) +
         7.83e-4 * f_GHz +
         5.88e-2;
}

G4ThreeVector G4CMPPhononBoundaryProcess::GetLambertianVector(G4ThreeVector surfaceNorm) {
  const G4int maxTries = 1000;
  G4int nTries = 0;
  do {
    reflectedKDir = G4CMP::LambertReflection(surfNorm);
  } while (nTries++ < maxTries &&
     !G4CMP::PhononVelocityIsInward(theLattice, mode,
            reflectedKDir, surfNorm));

  return reflectedKDir;
}

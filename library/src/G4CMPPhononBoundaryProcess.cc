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

#include "G4CMPPhononBoundaryProcess.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPSurfaceProperty.hh"
#include "G4CMPTrackInformation.hh"
#include "G4ExceptionSeverity.hh"
#include "G4GeometryTolerance.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4Navigator.hh"
#include "G4ParallelWorldProcess.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4TransportationManager.hh"
#include "G4VParticleChange.hh"
#include "G4LogicalSurface.hh"
#include "G4LogicalBorderSurface.hh"
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
  if (!IsGoodBoundary(aStep)) return &aParticleChange;

  if (verboseLevel>1) G4cout << GetProcessName() << "::PostStepDoIt" << G4endl;

  if (verboseLevel>2) {
    G4cout << " K direction: " << GetLocalWaveVector(aTrack).unit()
           << "\n P direction: " << aTrack.GetMomentumDirection() << G4endl;
  }

  ApplyBoundaryAction(aTrack, aStep, aParticleChange);
  return &aParticleChange;
}


G4bool G4CMPPhononBoundaryProcess::AbsorbTrack(const G4Track& aTrack,
                                               const G4Step& aStep) {
  G4double absMinK = matTable->GetConstProperty("absMinK");
  G4ThreeVector k = GetTrackInfo()->GetPhononK();
  
  return (G4CMPBoundaryUtils::AbsorbTrack(aTrack,aStep) &&
	  k*GetSurfaceNormal(aStep) > absMinK);
}


void
G4CMPPhononBoundaryProcess::DoReflection(const G4Track& aTrack, const G4Step& aStep,
					 G4ParticleChange& aParticleChange) {
  G4CMPTrackInformation* trackInfo = GetTrackInfo();

  if (verboseLevel>1) {
    G4cout << GetProcessName() << ": Track reflected "
           << trackInfo->GetReflectionCount() << " times." << G4endl;
  }

  // FIXME:  waveVector and reflectedKDir need to handle local/global rotations!
  G4ThreeVector waveVector = trackInfo->GetPhononK();
  G4int pol = GetPolarization(aStep.GetTrack());
  G4ThreeVector surfNorm = GetSurfaceNormal(aStep);

  if (verboseLevel>2)
    G4cout << " Old momentum direction " << waveVector.unit() << G4endl;

  G4double specProb = matTable->GetConstProperty("specProb");
  G4ThreeVector reflectedKDir;
  do {
    reflectedKDir = waveVector.unit();
    if (G4UniformRand() < specProb) {
      // Specular reflecton reverses momentum along normal
      G4double momNorm = reflectedKDir * surfNorm;
      reflectedKDir -= 2. * momNorm * surfNorm;
    } else {
      reflectedKDir = LambertReflection(surfNorm);
    }
  } while (!ReflectionIsGood(pol, reflectedKDir, surfNorm));

  if (verboseLevel>2)
    G4cout << " New momentum direction " << reflectedKDir << G4endl;

  G4ThreeVector vdir = theLattice->MapKtoVDir(pol, reflectedKDir);
  G4double v = theLattice->MapKtoV(pol, reflectedKDir);
  trackInfo->SetPhononK(reflectedKDir);
  aParticleChange.ProposeVelocity(v);
  aParticleChange.ProposeMomentumDirection(vdir);
}

G4bool G4CMPPhononBoundaryProcess::ReflectionIsGood(G4int polarization,
                                                    G4ThreeVector waveVector,
                                                    G4ThreeVector surfNorm) {
  G4ThreeVector vDir = theLattice->MapKtoVDir(polarization, waveVector);
  return vDir.dot(surfNorm) < 0.0;
}

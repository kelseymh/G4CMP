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
  : G4VPhononProcess(aName, fPhononReflection),
    kCarTolerance(G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()) {;}

G4double G4CMPPhononBoundaryProcess::
PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
                                     G4double previousStepSize,
                                     G4ForceCondition* condition) {
  return GetMeanFreePath(aTrack, previousStepSize, condition);
}

G4VParticleChange*
G4CMPPhononBoundaryProcess::PostStepDoIt(const G4Track& aTrack,
                                         const G4Step& aStep) {
  aParticleChange.Initialize(aTrack);
  G4StepPoint* postStepPoint = aStep.GetPostStepPoint();
  G4StepPoint* preStepPoint = aStep.GetPreStepPoint();

  // do nothing if the current step is not limited by a volume boundary,
  // or if it is the returning "null step" after a reflection
  if (postStepPoint->GetStepStatus() != fGeomBoundary ||
      aTrack.GetStepLength() <= kCarTolerance/2.) {
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }

  if (verboseLevel) {
    G4cout << GetProcessName() << "::PostStepDoIt length "
     << aTrack.GetStepLength() << G4endl;
  }

  G4VPhysicalVolume* thePrePV = preStepPoint->GetPhysicalVolume();
  G4VPhysicalVolume* thePostPV = postStepPoint->GetPhysicalVolume();

  if (thePrePV == thePostPV) {
    if (verboseLevel) {
      G4cerr << GetProcessName() << " ERROR: fGeomBoundary status set, but"
       << " pre- and post-step volumes are identical!" << G4endl;
    }
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }

  // do nothing if the current step is inbound from outside the original volume
  G4LatticePhysical* volLattice =
    G4LatticeManager::GetLatticeManager()->GetLattice(thePrePV);
  if (volLattice != theLattice) {
    if (verboseLevel>1) {
      G4cout << GetProcessName() << ": Track inbound after reflection"
       << G4endl;
    }
    return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);
  }

  if (verboseLevel>1) {
    G4cout <<   "    Track volume: " << aTrack.GetVolume()->GetName()
     << "\n  PreStep volume: " << thePrePV->GetName() << " @ "
     << preStepPoint->GetPosition()
     << "\n PostStep volume: " << thePostPV->GetName() << " @ "
     << postStepPoint->GetPosition()
     << G4endl;
  }

  // Grab surface information
  G4LogicalSurface* surfaceLog = G4LogicalBorderSurface::GetSurface(thePrePV,
                                                                    thePostPV);

  if (!surfaceLog) {
    G4Exception("G4CMPDriftBoundaryProcess::PostStepDoIt", "Boundary001",
                EventMustBeAborted, ("No surface defined between " +
                                    thePrePV->GetName() + " and " +
                                    thePostPV->GetName() + ".").c_str());
  }

  G4SurfaceProperty* surfProp = surfaceLog->GetSurfaceProperty();

  if (!surfProp) {
    G4Exception("G4CMPDriftBoundaryProcess::PostStepDoIt", "Boundary002",
                EventMustBeAborted, ("No surface property defined for " +
                                    surfaceLog->GetName() + ".").c_str());
  }

  if (verboseLevel>2) {
    G4cout << " K direction: " << GetLocalWaveVector(aTrack).unit()
           << "\n P direction: " << aTrack.GetMomentumDirection() << G4endl;
  }

  // If the particle doesn't get absorbed, it either reflects or transmits
  if (AbsorbTrack(aTrack, aStep, surfProp)) {
    return DoAbsorption(aTrack, aStep, surfProp);
  } else if (ReflectTrack(aTrack, aStep, surfProp)) {
    return DoReflection(aTrack, aStep, surfProp);
  } else {
    return DoTransmission(aTrack, aStep, surfProp);
  }

  return &aParticleChange;
}

void G4CMPPhononBoundaryProcess::LoadDataForTrack(const G4Track* track) {
  G4CMPProcessUtils::LoadDataForTrack(track);
  maxNumReflections = G4CMPConfigManager::GetMaxPhononBounces();
}


G4double G4CMPPhononBoundaryProcess::GetMeanFreePath(const G4Track&aTrack,
                                             G4double prevStepLength,
                                             G4ForceCondition* condition) {
  *condition = Forced;
  return DBL_MAX;
}

G4bool G4CMPPhononBoundaryProcess::AbsorbTrack(const G4Track& aTrack,
                                               const G4Step& aStep,
                                               const G4SurfaceProperty* surfProp) {
  // Check out this abomination:
  G4MaterialPropertiesTable*
    phonPropTable = const_cast<G4MaterialPropertiesTable*>(
                      static_cast<const G4CMPSurfaceProperty*>(surfProp)->
                        GetPhononMaterialPropertiesTablePointer()
                                                            );
  G4ThreeVector surfNorm = GetSurfaceNormal(aStep);
  G4double absProb = phonPropTable->GetConstProperty("absProb");
  G4double absMinK = phonPropTable->GetConstProperty("absMinK");

  G4ThreeVector k = GetTrackInfo()->GetPhononK();

  return (G4UniformRand() <= absProb && k*surfNorm > absMinK);
}

G4VParticleChange*
G4CMPPhononBoundaryProcess::DoAbsorption(const G4Track& aTrack,
                                         const G4Step& aStep,
                                         const G4SurfaceProperty* surfProp) {
  if (verboseLevel>1)
    G4cout << GetProcessName() << ": Track absorbed" << G4endl;

  G4double Ekin = aTrack.GetKineticEnergy();

  aParticleChange.ProposeNonIonizingEnergyDeposit(Ekin);
  aParticleChange.ProposeTrackStatus(fStopAndKill);
  return &aParticleChange;
}

G4bool G4CMPPhononBoundaryProcess::ReflectTrack(const G4Track& aTrack,
                                                const G4Step& aStep,
                                                const G4SurfaceProperty* surfProp) {
  G4MaterialPropertiesTable*
    phonPropTable = const_cast<G4MaterialPropertiesTable*>(
                        static_cast<const G4CMPSurfaceProperty*>(surfProp)->
                          GetPhononMaterialPropertiesTablePointer()
                                                            );
  G4double reflProb = phonPropTable->GetConstProperty("reflProb");
  return (G4UniformRand() <= reflProb);
}

G4VParticleChange*
G4CMPPhononBoundaryProcess::DoReflection(const G4Track& aTrack,
                                         const G4Step& aStep,
                                         const G4SurfaceProperty* surfProp) {
  const G4int maxRefl = G4CMPConfigManager::GetMaxPhononBounces();
  G4CMPTrackInformation* trackInfo = static_cast<G4CMPTrackInformation*>(
                            aTrack.GetAuxiliaryTrackInformation(fPhysicsModelID)
                                                                        );
  trackInfo->IncrementReflectionCount();
  if (trackInfo->GetReflectionCount() > maxRefl && maxRefl >= 0) {
    // If it reflects more than the user wants, we just kill it without
    // absorbing. This corresponds to the charge being too low energy to
    // detect.
    aParticleChange.ProposeTrackStatus(fStopAndKill);
    return &aParticleChange;
  }

  if (verboseLevel>1) {
    G4cout << GetProcessName() << ": Track reflected "
           << trackInfo->GetReflectionCount() << " times." << G4endl;
  }

  G4MaterialPropertiesTable*
    phonPropTable = const_cast<G4MaterialPropertiesTable*>(
                        static_cast<const G4CMPSurfaceProperty*>(surfProp)->
                          GetPhononMaterialPropertiesTablePointer()
                                                          );
  G4double specProb = phonPropTable->GetConstProperty("specProb");
  G4ThreeVector waveVector = trackInfo->GetPhononK();
  G4ThreeVector reflectedKDir = waveVector.unit();
  G4int pol = GetPolarization(aStep.GetTrack());
  G4ThreeVector surfNorm = GetSurfaceNormal(aStep);

  if (verboseLevel>2)
    G4cout << " Old momentum direction " << waveVector.unit() << G4endl;

  do {
    if (G4UniformRand() < specProb) {
      // Specular reflecton reverses momentum along normal
      G4double momNorm = reflectedKDir * surfNorm;
      reflectedKDir -= 2. * momNorm * surfNorm;
    } else {
        reflectedKDir = LambertReflection(surfNorm);
    }
    waveVector = waveVector.mag()*reflectedKDir;
  } while (!ReflectionIsGood(pol, waveVector, surfNorm));

  if (verboseLevel>2)
    G4cout << " New momentum direction " << reflectedKDir << G4endl;

  G4ThreeVector vdir = theLattice->MapKtoVDir(pol, waveVector);
  G4double v = theLattice->MapKtoV(pol, waveVector);
  trackInfo->SetPhononK(waveVector);
  aParticleChange.ProposeVelocity(v);
  aParticleChange.ProposeMomentumDirection(vdir);

  return &aParticleChange;
}

G4VParticleChange*
G4CMPPhononBoundaryProcess::DoTransmission(const G4Track& aTrack,
                                           const G4Step& aStep,
                                           const G4SurfaceProperty* surfProp) {
  //noop
  return &aParticleChange;
}

G4bool G4CMPPhononBoundaryProcess::ReflectionIsGood(G4int polarization,
                                                    G4ThreeVector waveVector,
                                                    G4ThreeVector surfNorm) {
  G4ThreeVector vDir = theLattice->MapKtoVDir(polarization, waveVector);
  return vDir.dot(surfNorm) < 0.0;
}

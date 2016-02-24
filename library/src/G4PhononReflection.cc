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

#include "G4CMPConfigManager.hh"
#include "G4CMPSurfaceProperty.hh"
#include "G4CMPTrackInformation.hh"
#include "G4PhononReflection.hh"
#include "G4ExceptionSeverity.hh"
#include "G4GeometryTolerance.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4PhononTrackMap.hh"
#include "G4PhononPolarization.hh"
#include "G4PhononLong.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononTransSlow.hh"
#include "G4PhysicalConstants.hh"
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


G4PhononReflection::G4PhononReflection(const G4String& aName)
  : G4VPhononProcess(aName, fPhononReflection),
    kCarTolerance(G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()) {;}

G4VParticleChange* G4PhononReflection::PostStepDoIt(const G4Track& aTrack,
                                                    const G4Step& aStep) {
  aParticleChange.Initialize(aTrack);
  G4StepPoint* postStepPoint = aStep.GetPostStepPoint();
  G4StepPoint* preStepPoint = aStep.GetPreStepPoint();

  // Do nothing but return if the step is not boundary limited, or if
  // the step is too short. This is to allow actual reflection where after
  // the first boundary crossing a second, infinitesimal, step occurs crossing
  // back into the original volume
  if (postStepPoint->GetStepStatus() != fGeomBoundary ||
      aTrack.GetStepLength() <= kCarTolerance/2.) {
    return &aParticleChange;
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
  // FIXME: Isn't this test redundant with the kCarTolerance test?
  G4LatticePhysical*
    volLattice = G4LatticeManager::GetLatticeManager()->GetLattice(thePrePV);
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
  G4LogicalSurface* surface = G4LogicalBorderSurface::GetSurface(thePrePV,
                                                                 thePostPV);
  G4CMPSurfaceProperty* borderSurface = nullptr;
  if (surface) {
    borderSurface =
      static_cast<G4CMPSurfaceProperty*>(surface->GetSurfaceProperty());
  } else {
    if (verboseLevel>1) {
      G4cerr << GetProcessName() << ": No border surface defined for "
             << thePrePV->GetName() << " to "  << thePostPV->GetName()
             << G4endl;
    }

    return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);
  }

  if (!borderSurface) {
    if (verboseLevel>1) {
      G4cerr << GetProcessName() << ": No surface properties defined for "
             << thePrePV->GetName() << " to "  << thePostPV->GetName()
             << G4endl;
    }

    return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);
  }

  // Grab phonon properties from the surface
  G4MaterialPropertiesTable*
    phonPropTable = const_cast<G4MaterialPropertiesTable*>(
      borderSurface->GetPhononMaterialPropertiesTablePointer());
  absProb = phonPropTable->GetConstProperty("absProb");
  reflProb = phonPropTable->GetConstProperty("reflProb");
  specProb = phonPropTable->GetConstProperty("specProb");
  absMinK = phonPropTable->GetConstProperty("minK");


  // Get outward normal using G4Navigator method (more reliable than G4VSolid)
  G4int navID = G4ParallelWorldProcess::GetHypNavigatorID();
  std::vector<G4Navigator*>::iterator iNav =
    G4TransportationManager::GetTransportationManager()->GetActiveNavigatorsIterator();
  G4bool goodNorm;
  surfNorm =
    iNav[navID]->GetGlobalExitNormal(postStepPoint->GetPosition(), &goodNorm);
  if (!goodNorm) {
    G4cerr << GetProcessName() << " ERROR:  Cannot get normal at surface of "
           << thePrePV->GetName() << " @ " << postStepPoint->GetPosition()
           << G4endl;
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  } else if (verboseLevel>2) {
    G4cout << " Normal " << surfNorm << " @ " << postStepPoint->GetPosition()
           << G4endl;
  }

  trackInfo = static_cast<G4CMPTrackInformation*>(
                          aTrack.GetAuxiliaryTrackInformation(fPhysicsModelID));
  waveVector = trackInfo->GetPhononK();


  // If the particle doesn't get absorbed, it either reflects or transmits
  if (AbsorbTrack(aStep)) {
    return DoAbsorption(aStep);
  } else if (ReflectTrack(aStep)) {
    return DoReflection(aStep);
  } // else the particle passes through the boundary

  return &aParticleChange; 
}

G4double G4PhononReflection::GetMeanFreePath(const G4Track&aTrack,
                                             G4double prevStepLength,
                                             G4ForceCondition* condition) {
  // Always return DBL_MAX and Forced. This ensures that the process is
  // called at the end of every step. In PostStepDoIt the process
  // decides whether the step encountered a volume boundary and a
  // reflection should be applied
  *condition = Forced;
  return DBL_MAX;
}

G4bool G4PhononReflection::AbsorbTrack(const G4Step& aStep) {
  return (G4UniformRand() <= absProb && waveVector*surfNorm > absMinK);
}

G4VParticleChange* G4PhononReflection::DoAbsorption(const G4Step& aStep) {
  if (verboseLevel>1)
    G4cout << GetProcessName() << ": Track absorbed" << G4endl;

  G4Track* aTrack = aStep.GetTrack();
  G4double Ekin = aTrack->GetKineticEnergy();

  aParticleChange.ProposeNonIonizingEnergyDeposit(Ekin);
  aParticleChange.ProposeTrackStatus(fStopAndKill);
  return &aParticleChange;
}

G4bool G4PhononReflection::ReflectTrack(const G4Step& aStep) {
  return (G4UniformRand() <= reflProb);
}

G4VParticleChange* G4PhononReflection::DoReflection(const G4Step& aStep) {
  const G4int maxRefl = G4CMPConfigManager::GetMaxPhononBounces();
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

  if (verboseLevel>2)
    G4cout << " Old momentum direction " << waveVector.unit() << G4endl;


  G4ThreeVector reflectedKDir = waveVector.unit();
  G4int pol = GetPolarization(aStep.GetTrack());
  do {
    if (G4UniformRand() < specProb) {
      // Specular reflecton reverses momentum along normal
      G4double momNorm = reflectedKDir * surfNorm;
      reflectedKDir -= 2. * momNorm * surfNorm;
    } else {
        reflectedKDir = LambertReflection(surfNorm);
    }
    waveVector = waveVector.mag()*reflectedKDir;
  } while (!ReflectionIsGood(pol));

  if (verboseLevel>2)
    G4cout << " New momentum direction " << reflectedKDir << G4endl;

  G4ThreeVector vdir = theLattice->MapKtoVDir(pol, waveVector);
  G4double v = theLattice->MapKtoV(pol, waveVector);
  trackInfo->SetPhononK(waveVector);
  aParticleChange.ProposeVelocity(v);
  aParticleChange.ProposeMomentumDirection(vdir);

  return &aParticleChange;
}

G4bool G4PhononReflection::ReflectionIsGood(G4int polarization) {
  G4ThreeVector vDir = theLattice->MapKtoVDir(polarization, waveVector);
  return vDir.dot(surfNorm) < 0.0;
}

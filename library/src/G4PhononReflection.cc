//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
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
  // This process handles the interaction of phonons with
  // boundaries. Implementation of this class is highly geometry
  // dependent.Currently, phonons are killed when they reach a
  // boundary. If the other side of the boundary was Al, a hit is
  // registered.
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
  // FIXME: Shouldn't this test be redundant with the kCarTolerance test?
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

  // Get outward normal using G4Navigator method (more reliable than G4VSolid)
  G4int navID = G4ParallelWorldProcess::GetHypNavigatorID();
  std::vector<G4Navigator*>::iterator iNav =
    G4TransportationManager::GetTransportationManager()->GetActiveNavigatorsIterator();

  G4bool goodNorm;
  G4ThreeVector surfNorm =
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

  const G4int maxRefl = G4CMPConfigManager::GetMaxPhononBounces();
  G4CMPTrackInformation* trackInfo = static_cast<G4CMPTrackInformation*>(
    aTrack.GetAuxiliaryTrackInformation(fPhysicsModelID));
  G4int nRefl = trackInfo->GetReflectionCount();
  if (maxRefl<0 || nRefl < maxRefl) {
    trackInfo->IncrementReflectionCount();
    return DoReflection(aStep, surfNorm, borderSurface, trackInfo);
  } else {
    if (verboseLevel) {
      G4cout << GetProcessName() << " WARNING: Phonon has reflected "
             << maxRefl << " times. Track being killed." << G4endl;
    }
    G4double eKin = aTrack.GetKineticEnergy();
    aParticleChange.ProposeNonIonizingEnergyDeposit(eKin);
    aParticleChange.ProposeTrackStatus(fStopAndKill);
  }

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

G4VParticleChange*
G4PhononReflection::DoReflection(const G4Step& aStep,
                                 const G4ThreeVector& surfNorm,
                                 const G4CMPSurfaceProperty* surfProp,
                                 G4CMPTrackInformation* trackInfo) {
  if (verboseLevel>1)
    G4cout << GetProcessName() << ": Track reflected" << G4endl;

  G4ThreeVector k = trackInfo->GetPhononK();
  if (verboseLevel>2)
    G4cout << " Old momentum direction " << k.unit() << G4endl;

  // Grab phonon properties from the surface
  G4MaterialPropertiesTable*
    phonPropTable = const_cast<G4MaterialPropertiesTable*>(
      surfProp->GetPhononMaterialPropertiesTablePointer());
  G4double absProb       = phonPropTable->GetConstProperty("absProb");
  G4double specProb      = phonPropTable->GetConstProperty("specProb");

  if (G4UniformRand() < absProb) {
    // Kill the track, but potentially make new reflected phonons.
    aParticleChange.ProposeTrackStatus(fStopAndKill);

    G4double Ekin = aStep.GetPostStepPoint()->GetKineticEnergy();
    // Track deposits some energy and may spawn new phonons.
    std::vector<G4double> phononEnergies;
    // NOTE: phononEnergies mutates in KaplanPhononQP
    G4double EDep = KaplanPhononQP(Ekin, phonPropTable, phononEnergies);
    aParticleChange.ProposeNonIonizingEnergyDeposit(EDep);

    for (G4double E : phononEnergies) {
      // TODO: Map E to K properly.
      G4double kmag = k.mag()*E/Ekin;
      G4ThreeVector reflectedKDir = k.unit();
      G4int pol = ChoosePolarization();
      do {
        reflectedKDir = LambertReflection(surfNorm);
      } while (!ReflectionIsGood(pol, kmag*reflectedKDir, surfNorm));

      CreatePhonon(pol,kmag*reflectedKDir,E);
    }
  } else {
    G4ThreeVector reflectedKDir = k.unit();
    G4int pol = GetPolarization(aStep.GetTrack());
    do {
      if (G4UniformRand() < specProb) {
        // Specular reflecton reverses momentum along normal
        G4double momNorm = reflectedKDir * surfNorm;
        reflectedKDir -= 2. * momNorm * surfNorm;
      } else {
        reflectedKDir = LambertReflection(surfNorm);
      }
      k = k.mag()*reflectedKDir;
    } while (!ReflectionIsGood(pol, k, surfNorm));

    if (verboseLevel>2)
      G4cout << " New momentum direction " << reflectedKDir << G4endl;

    G4ThreeVector vdir = theLattice->MapKtoVDir(pol, k);
    G4double v = theLattice->MapKtoV(pol, k);
    trackInfo->SetPhononK(k);
    aParticleChange.ProposeVelocity(v);
    aParticleChange.ProposeMomentumDirection(vdir);
  } // if absorbed

  return &aParticleChange;
}

G4bool G4PhononReflection::ReflectionIsGood(G4int polarization,
                                            const G4ThreeVector& k,
                                            const G4ThreeVector& surfNorm) {
  G4ThreeVector vDir = theLattice->MapKtoVDir(polarization, k);
  return vDir.dot(surfNorm) < 0.0;
}

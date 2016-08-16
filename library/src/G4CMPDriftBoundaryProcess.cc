/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
//
// 20140331  Inherit from G4CMPVDriftProcess to get subtype enforcement
// 20141029  Get output hits file from configuration manager
// 20150122  Use verboseLevel instead of compiler flag for debugging
// 20150212  Remove file IO. Use sensitive detectors instead
// 20150603  Add functionality to globally limit reflections

#include "G4CMPDriftBoundaryProcess.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4CMPSurfaceProperty.hh"
#include "G4GeometryTolerance.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4ParallelWorldProcess.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TransportationManager.hh"
#include "G4VParticleChange.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"
#include <vector>


G4CMPDriftBoundaryProcess::G4CMPDriftBoundaryProcess(const G4String& name,
                                                     G4CMPProcessSubType type)
  : G4CMPVDriftProcess(name, type),
    kCarTolerance(G4GeometryTolerance::GetInstance()->GetSurfaceTolerance())
{}

G4double G4CMPDriftBoundaryProcess::
PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
				     G4double previousStepSize,
				     G4ForceCondition* condition) {
  return GetMeanFreePath(aTrack, previousStepSize, condition);
}

G4double G4CMPDriftBoundaryProcess::
GetMeanFreePath(const G4Track& /*aTrack*/,G4double /*previousStepSize*/,
		G4ForceCondition* condition) {
  *condition = Forced;
  return DBL_MAX;
}

G4VParticleChange* 
G4CMPDriftBoundaryProcess::PostStepDoIt(const G4Track& aTrack,
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
    if (aTrack.GetDefinition() == G4CMPDriftElectron::Definition()) {
      G4cout << " K_valley (" << GetValleyIndex(aTrack) << ") direction: "
	     << theLattice->MapPtoK_valley(GetValleyIndex(aTrack),
					   GetLocalMomentum(aTrack)).unit()
	     << G4endl;
    }
    G4cout << " K direction: " << GetLocalWaveVector(aTrack).unit()
           << "\n P direction: " << GetLocalMomentum(aTrack).unit() << G4endl;
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

void G4CMPDriftBoundaryProcess::LoadDataForTrack(const G4Track* track) {
  G4CMPVDriftProcess::LoadDataForTrack(track);
  maxNumReflections = G4CMPConfigManager::GetMaxChargeBounces();
  numReflections = 0;		// New track, no bounces
}

// Decide and apply different surface actions; subclasses may override

G4bool G4CMPDriftBoundaryProcess::AbsorbTrack(const G4Track& aTrack,
                                              const G4Step& aStep,
                                              const G4SurfaceProperty* surfProp)
                                              const {
  // Check out this abomination:
  G4MaterialPropertiesTable*
    chargePropTable = const_cast<G4MaterialPropertiesTable*>(
                        static_cast<const G4CMPSurfaceProperty*>(surfProp)->
                          GetChargeMaterialPropertiesTablePointer()
                                                            );
  G4double absProb = chargePropTable->GetConstProperty("absProb");
  G4double absMinK;
  if (aTrack.GetDefinition() == G4CMPDriftElectron::Definition()) {
    absMinK = chargePropTable->GetConstProperty("minKElec");
  } else if (aTrack.GetDefinition() == G4CMPDriftHole::Definition()) {
    absMinK = chargePropTable->GetConstProperty("minKHole");
  } else {
    G4Exception("G4CMPDriftBoundaryProcess::AbsorbTrack", "Boundary003",
                EventMustBeAborted, "Invalid particle for this process.");
  }

  G4ThreeVector kvec = GetLocalWaveVector(aTrack);

  // NOTE:  K vector above is in local coords, must use local normal
  G4ThreeVector surfNorm = GetLocalDirection(GetSurfaceNormal(aStep));

  if (verboseLevel>2) {
    G4cout << " AbsorbTrack: absProb " << absProb
	   << " local k-perp " << kvec*surfNorm
	   <<" >? absMinK " << absMinK << G4endl;
  }

  return (G4UniformRand() <= absProb && kvec*surfNorm > absMinK);
}

G4VParticleChange*
G4CMPDriftBoundaryProcess::DoAbsorption(const G4Track& aTrack,
                                        const G4Step& /*aStep*/,
                                        const G4SurfaceProperty* /*surfProp*/) {
  if (verboseLevel>1)
    G4cout << GetProcessName() << ": Track absorbed" << G4endl;

  G4double Ekin = GetKineticEnergy(aTrack);

  aParticleChange.ProposeNonIonizingEnergyDeposit(Ekin);
  aParticleChange.ProposeTrackStatus(fStopAndKill);
  return &aParticleChange;
}

G4bool G4CMPDriftBoundaryProcess::ReflectTrack(const G4Track& /*aTrack*/,
                                               const G4Step& /*aStep*/,
                                               const G4SurfaceProperty* surfProp)
                                               const {
  G4MaterialPropertiesTable*
    chargePropTable = const_cast<G4MaterialPropertiesTable*>(
                        static_cast<const G4CMPSurfaceProperty*>(surfProp)->
                          GetChargeMaterialPropertiesTablePointer()
                                                            );
  G4double reflProb = chargePropTable->GetConstProperty("reflProb");

  if (verboseLevel>2)
    G4cout << " ReflectTrack: reflProb " << reflProb << G4endl;

  return (G4UniformRand() <= reflProb);
}

G4VParticleChange*  
G4CMPDriftBoundaryProcess::DoReflection(const G4Track& aTrack,
                                        const G4Step& aStep,
                                        const G4SurfaceProperty* /*surfProp*/) {
  if (++numReflections > maxNumReflections && maxNumReflections >= 0) {
    // If it reflects more than the user wants, we just kill it without
    // absorbing.
    if (verboseLevel>1)
      G4cout << GetProcessName() << ": Track reflected more than "
             << maxNumReflections << " times. Track killed." << G4endl;

    aParticleChange.ProposeTrackStatus(fStopAndKill);
    return &aParticleChange;
  }

  if (verboseLevel>1)
    G4cout << GetProcessName() << ": Track reflected" << G4endl;

  G4ThreeVector surfNorm = GetSurfaceNormal(aStep);

  // Electrons and holes need to be handled separately until we further
  // generalize the physics.

  if (aTrack.GetDefinition() == G4CMPDriftElectron::Definition()) {
    G4ThreeVector vel = GetGlobalVelocityVector(aTrack);

    if (verboseLevel>2)
      G4cout << " Old velocity direction " << vel.unit() << G4endl;

    // Specular reflecton reverses velocity along normal
    G4double velNorm = vel * surfNorm;
    vel -= 2.*velNorm*surfNorm;

    if (verboseLevel>2)
      G4cout << " New velocity direction " << vel.unit() << G4endl;

    // Convert velocity back to momentum and update direction
    RotateToLocalDirection(vel);
    G4ThreeVector p = theLattice->MapV_elToP(GetCurrentValley(), vel);
    RotateToGlobalDirection(p);

    if (verboseLevel>2) {
      G4cout << " New momentum direction " << p.unit() << G4endl;

      // SANITY CHECK:  Does new momentum get back to new velocity?
      G4ThreeVector vnew = theLattice->MapPtoV_el(GetCurrentValley(),
                                                  GetLocalDirection(p));
      RotateToGlobalDirection(vnew);
      G4cout << " Cross-check new v dir  " << vnew.unit() << G4endl;
    }

    FillParticleChange(GetCurrentValley(), p);	// Handle effective mass, vel
  } else if (aTrack.GetDefinition() == G4CMPDriftHole::Definition()) {
    G4ThreeVector momDir = aStep.GetPostStepPoint()->GetMomentumDirection();
    if (verboseLevel>2)
      G4cout << " Old momentum direction " << momDir << G4endl;

    G4double momNorm = momDir * surfNorm;
    momDir -= 2.*momNorm*surfNorm;

    if (verboseLevel>2)
      G4cout << " New momentum direction " << momDir << G4endl;

    aParticleChange.ProposeMomentumDirection(momDir);
  } else {
    G4Exception("G4CMPDriftBoundaryProcess::DoReflection", "Boundary004",
                EventMustBeAborted, "Invalid particle for this process.");
  }

  return &aParticleChange;
}

G4VParticleChange*
G4CMPDriftBoundaryProcess::DoTransmission(const G4Track& /*aTrack*/,
                                          const G4Step& /*aStep*/,
                                          const G4SurfaceProperty* /*surfProp*/) {
  if (verboseLevel>1)
    G4cout << GetProcessName() << ": Track transmitted" << G4endl;

  //noop - "Move along, Particle."
  return &aParticleChange;
}

G4ThreeVector
G4CMPDriftBoundaryProcess::GetSurfaceNormal(const G4Step& aStep) const {
  // Get outward normal using G4Navigator method (more reliable than G4VSolid)
  G4int navID = G4ParallelWorldProcess::GetHypNavigatorID();
  std::vector<G4Navigator*>::iterator iNav =
    G4TransportationManager::GetTransportationManager()->GetActiveNavigatorsIterator();

  G4bool goodNorm;
  G4ThreeVector surfNorm = iNav[navID]->GetGlobalExitNormal(
                                      aStep.GetPostStepPoint()->GetPosition(),
                                      &goodNorm);
  if (!goodNorm) {
    G4VPhysicalVolume* thePrePV = aStep.GetPreStepPoint()->GetPhysicalVolume();
    G4VPhysicalVolume* thePostPV = aStep.GetPostStepPoint()->GetPhysicalVolume();
    G4Exception("G4CMPDriftBoundaryProcess::PostStepDoIt", "Boundary001",
                EventMustBeAborted, ("Can't get normal vector of surface between " +
                                    thePrePV->GetName() + " and " +
                                    thePostPV->GetName()+ ".").c_str());
  } else if (verboseLevel>2) {
    G4cout << " Normal " << surfNorm << " @ "
           << aStep.GetPostStepPoint()->GetPosition() << G4endl;
  }

  return surfNorm;
}

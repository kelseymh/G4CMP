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
#include "G4CMPMeshElectricField.hh"
#include "G4CMPSurfaceProperty.hh"
#include "G4CMPTrackInformation.hh"
#include "G4Field.hh"
#include "G4FieldManager.hh"
#include "G4GeometryTolerance.hh"
#include "G4LatticeLogical.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4ParallelWorldProcess.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhononPolarization.hh"
#include "G4RandomDirection.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4VParticleChange.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"
#include <vector>


G4CMPDriftBoundaryProcess::G4CMPDriftBoundaryProcess(const G4String& name)
  : G4CMPVDriftProcess("G4CMPDriftBoundary", fChargeBoundary),
    kCarTolerance(G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()),
    shortName(name), theCarrier(nullptr), reflProb(0.), absProb(1.),
    absMinK(0.), numberOfReflections(0)
{;}

G4CMPDriftBoundaryProcess::~G4CMPDriftBoundaryProcess() {;}


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
  G4LogicalSurface* surface = G4LogicalBorderSurface::GetSurface(thePrePV,
                                                                 thePostPV);
  G4CMPSurfaceProperty* borderSurface;
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

  G4MaterialPropertiesTable*
  chargePropTable = const_cast<G4MaterialPropertiesTable*>(
                      borderSurface->GetChargeMaterialPropertiesTablePointer());
  reflProb = chargePropTable->GetConstProperty("reflProb");
  absProb = chargePropTable->GetConstProperty("absProb");

  if (theCarrier == G4CMPDriftElectron::Definition()) {
    absMinK = chargePropTable->GetConstProperty("minKElec");
  } else if (theCarrier == G4CMPDriftHole::Definition()) {
    absMinK = chargePropTable->GetConstProperty("minKHole");
  } else {
    G4Exception("G4CMPDriftBoundaryProcess::PostStepDoIt", "Boundary001",
                EventMustBeAborted, "Invalid particle for this process.");
  }

  // Get outward normal using G4Navigator method (more reliable than G4VSolid)
  G4int navID = G4ParallelWorldProcess::GetHypNavigatorID();
  std::vector<G4Navigator*>::iterator iNav =
    G4TransportationManager::GetTransportationManager()->GetActiveNavigatorsIterator();

  G4bool goodNorm;
  surfNorm = iNav[navID]->GetGlobalExitNormal(postStepPoint->GetPosition(),
                                              &goodNorm);
  if (!goodNorm) {
    G4cerr << GetProcessName() << " ERROR:  Cannot get normal at surface of "
           << thePrePV->GetName() << " @ " << postStepPoint->GetPosition()
           << G4endl;
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  } else if (verboseLevel>2) {
    G4cout << " Normal " << surfNorm << " @ " << postStepPoint->GetPosition()
           << G4endl;
  }

  // If the particle doesn't get absorbed, it either reflects or transmits
  if (AbsorbTrack(aStep)) {
    return DoAbsorption(aStep);
  } else if (ReflectTrack(aStep)) {
    return DoReflection(aStep);
  } // else the particle passes through the boundary

  if (verboseLevel>2) {
    G4cout << " K direction: " << GetLocalWaveVector(aTrack).unit()
           << "\n P direction: " << aTrack.GetMomentumDirection() << G4endl;
  }

  return &aParticleChange;
}

void G4CMPDriftBoundaryProcess::LoadDataForTrack(const G4Track* track) {
  G4CMPVDriftProcess::LoadDataForTrack(track);
  theCarrier = track->GetDefinition();
  maxRefl = G4CMPConfigManager::GetMaxChargeBounces();
  numberOfReflections = 0;		// New track, no bounces
}

// Decide and apply different surface actions; subclasses may override

G4bool G4CMPDriftBoundaryProcess::AbsorbTrack(const G4Step& aStep) {
  // Universal absorption threshold
  return (G4UniformRand() <= absProb &&
          GetLocalWaveVector(aStep.GetTrack())*surfNorm > absMinK);
}

G4VParticleChange*
G4CMPDriftBoundaryProcess::DoAbsorption(const G4Step& aStep) {
  if (verboseLevel>1)
    G4cout << GetProcessName() << ": Track absorbed" << G4endl;

  G4Track* aTrack = aStep.GetTrack();
  G4double Ekin = GetKineticEnergy(aTrack);

  aParticleChange.ProposeNonIonizingEnergyDeposit(Ekin);
  aParticleChange.ProposeTrackStatus(fStopAndKill);
  return &aParticleChange;
}

G4bool G4CMPDriftBoundaryProcess::ReflectTrack(const G4Step& aStep) {
  if (verboseLevel>1) {
    G4cout << GetProcessName() << ": Track reflected " << numberOfReflections
           << " times." << G4endl;
  }
  return (G4UniformRand() <= reflProb);
}

G4VParticleChange*  
G4CMPDriftBoundaryProcess::DoReflection(const G4Step& aStep) {
  ++numberOfReflections;
  if (numberOfReflections > maxRefl && maxRefl >= 0) {
    // If it reflects more than the user wants, we just kill it without
    // absorbing. This corresponds to the charge being too low energy to
    // detect.
    aParticleChange.ProposeTrackStatus(fStopAndKill);
    return &aParticleChange;
  }

  if (verboseLevel>1)
    G4cout << GetProcessName() << ": Track reflected" << G4endl;

  G4ThreeVector momDir = aStep.GetPostStepPoint()->GetMomentumDirection();
  if (verboseLevel>2)
    G4cout << " Old momentum direction " << momDir << G4endl;

  // Electrons and holes need to be handled separately until we further
  // generalize the physics.

  if (theCarrier == G4CMPDriftElectron::Definition()) {
    G4Track* aTrack = aStep.GetTrack();
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
  } else if (theCarrier == G4CMPDriftHole::Definition()) {
    G4double momNorm = momDir * surfNorm;
    momDir -= 2.*momNorm*surfNorm;

    if (verboseLevel>2)
      G4cout << " New momentum direction " << momDir << G4endl;

    aParticleChange.ProposeMomentumDirection(momDir);
  } else {
    G4Exception("G4CMPDriftBoundaryProcess::DoReflection", "Boundary002",
                EventMustBeAborted, "Invalid particle for this process.");
  }

  return &aParticleChange;
}

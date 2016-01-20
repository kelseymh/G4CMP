// $Id$
//
// 20140331  Inherit from G4CMPVDriftProcess to get subtype enforcement
// 20141029  Get output hits file from configuration manager
// 20150122  Use verboseLevel instead of compiler flag for debugging
// 20150212  Remove file IO. Use sensitive detectors instead
// 20150603  Add functionality to globally limit reflections

#include "G4CMPVDriftBoundaryProcess.hh"
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


G4CMPVDriftBoundaryProcess::G4CMPVDriftBoundaryProcess(const G4String& name,
                                         const G4ParticleDefinition* carrier)
  : G4CMPVDriftProcess("G4CMP"+name+"Boundary", fChargeBoundary),
    kCarTolerance(G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()),
    theCarrier(carrier), shortName(name), absProb(1.), absDeltaV(0.), 
    absMinKElec(0.), absMinKHole(0.), electrodeV(0.), numberOfReflections(0)
{;}

G4CMPVDriftBoundaryProcess::~G4CMPVDriftBoundaryProcess() {;}


G4bool 
G4CMPVDriftBoundaryProcess::IsApplicable(const G4ParticleDefinition& aPD) {
  return (&aPD==theCarrier);
}


G4double G4CMPVDriftBoundaryProcess::
PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
				     G4double previousStepSize,
				     G4ForceCondition* condition) {
  return GetMeanFreePath(aTrack, previousStepSize, condition);
}

G4double G4CMPVDriftBoundaryProcess::
GetMeanFreePath(const G4Track& /*aTrack*/,G4double /*previousStepSize*/,
		G4ForceCondition* condition) {
  *condition = Forced;
  return DBL_MAX;
}

G4VParticleChange* 
G4CMPVDriftBoundaryProcess::PostStepDoIt(const G4Track& aTrack,
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

  G4MaterialPropertiesTable* chargePropTable =
    const_cast<G4MaterialPropertiesTable*>(
      borderSurface->GetChargeMaterialPropertiesTablePointer());
  absProb = chargePropTable->GetConstProperty("absProb");
  electrodeV = chargePropTable->GetConstProperty("electrodeV");
  absDeltaV = chargePropTable->GetConstProperty("absDeltaV");
  absMinKElec = chargePropTable->GetConstProperty("minKElec");
  absMinKHole = chargePropTable->GetConstProperty("minKHole");

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

  // Test #1 and #2: Abosrption at surface without signal
  if (AbsorbTrack(aStep)) {
    return DoAbsorption(aStep);
  }

  if (verboseLevel>2) {
    G4cout <<   " K direction: " << GetLocalWaveVector(aTrack).unit()
	   << "\n P direction: " << aTrack.GetMomentumDirection() << G4endl;
  }

  // Test #3: If landed on an electrode.
  if (HitElectrode(aStep)) {
    return DoElectrodeHit(aStep);
  }

  // Test #4: Charge carrier reflects off surface
  const G4int maxRefl = G4CMPConfigManager::GetMaxChargeBounces();
  if (maxRefl<0 || numberOfReflections < maxRefl) {
    if (ReflectTrack(aStep)) {
      numberOfReflections++;
      return DoReflection(aStep);
    }
  } else {
    if (verboseLevel>1) {
      G4cout << GetProcessName() << ": Track reflected " << numberOfReflections
	     << " times." << G4endl;
    }
    return DoAbsorption(aStep);
  }

  // Should never reach this code!
  if (verboseLevel>1) {
    G4cerr << GetProcessName() << " ERROR: Track boundary process failed"
	   << G4endl;
  } 

  return &aParticleChange;
}

void G4CMPVDriftBoundaryProcess::LoadDataForTrack(const G4Track* track) {
  if (track->GetDefinition() != theCarrier) {
    G4cerr << GetProcessName() << " ERROR:  Track type "
     << track->GetDefinition()->GetParticleName() << " not valid" << G4endl;
    return;
  }

  G4CMPVDriftProcess::LoadDataForTrack(track);
  numberOfReflections = 0;		// New track, no bounces
}

// Decide and apply different surface actions; subclasses may override

G4bool G4CMPVDriftBoundaryProcess::AbsorbTrack(const G4Step& aStep) {
  // Universal absorption threshold
  return (G4UniformRand() <= absProb);
}

G4VParticleChange* 
G4CMPVDriftBoundaryProcess::DoAbsorption(const G4Step& aStep) {
  if (verboseLevel>1)
    G4cout << GetProcessName() << ": Track absorbed" << G4endl;

  G4Track* aTrack = aStep.GetTrack();
  G4double Ekin = GetKineticEnergy(aTrack);

/* FIXME: Eventually uncomment this code to produce a real phonon that will
   propagate into the next volume.

  // Create real phonon to be propagated, with random polarization
  static const G4double genPhonon = G4CMPConfigManager::GetGenPhonons();
  if (genPhonon > 0. && G4UniformRand() < genPhonon) {
    G4int polarization = ChoosePolarization();

    // Conserve energy and momentum(mag)
    // As a first guess, let's just take the carriers k and give it to the
    // phonon. We have to make sure if gives us a v that points out of the volume.
    G4ThreeVector k = GetLocalWaveVector(aTrack);
    G4double kmag = k.mag();
    G4ThreeVector vdir = theLattice->MapKtoVDir(polarization, k);
    RotateToGlobalDirection(vdir);
    // FIXME: We should have a MapVtoK map to make this not slow.
    while (vdir*surfNorm < 0.5) {
      G4ThreeVector kdir = G4RandomDirection();
      k = kmag*kdir;
      MakeLocalPhononK(k);
      vdir = theLattice->MapKtoVDir(polarization, k);
      RotateToGlobalDirection(vdir);
    }

    G4Track* phonon = CreatePhonon(polarization, k, Ekin);
    aParticleChange.SetNumberOfSecondaries(1);
    aParticleChange.AddSecondary(phonon);
  }
*/

  aParticleChange.ProposeNonIonizingEnergyDeposit(Ekin);
  aParticleChange.ProposeTrackStatus(fStopAndKill);
  return &aParticleChange;
}

G4bool G4CMPVDriftBoundaryProcess::HitElectrode(const G4Step& aStep) {
  G4VPhysicalVolume* thePrePV = aStep.GetPreStepPoint()->GetPhysicalVolume();
  G4FieldManager* fMan = thePrePV->GetLogicalVolume()->GetFieldManager();

  if (!fMan || !fMan->DoesFieldExist()) {	// Sanity check, must have field
    G4cerr << GetProcessName() << " ERROR: No Field for volume "
	   << thePrePV->GetName() << G4endl;
    return false;
  }

  // Can only get potential (to identify electrode) from G4CMP special field
  const G4CMPMeshElectricField* field = 
    dynamic_cast<const G4CMPMeshElectricField*>(fMan->GetDetectorField());
  if (!field) return false;

  // Extract potential at hit point from G4CMP-specific field interface
  G4double posVec[4] = { 4*0. };
  GetLocalPosition(aStep.GetTrack(), posVec);
  G4double potential = field->GetPotential(posVec);
  
  return (fabs(potential-electrodeV) <= absDeltaV);
}

// Default behaviour just kills the track, may need to have "state flag"
// akin to G4OpBoundaryProcess.

G4VParticleChange* 
G4CMPVDriftBoundaryProcess::DoElectrodeHit(const G4Step& aStep) {
  if (verboseLevel>1)
    G4cout << GetProcessName() << ": Track hit electrode" << G4endl;

  return DoAbsorption(aStep);
}


// Default behaviour assumes normal, scalar kinematics

G4bool G4CMPVDriftBoundaryProcess::ReflectTrack(const G4Step& aStep) {
  // Track is outbound hitting surface of volume
  return (aStep.GetPostStepPoint()->GetMomentumDirection()*surfNorm > 0.);
}

G4VParticleChange*  
G4CMPVDriftBoundaryProcess::DoReflection(const G4Step& aStep) {
  if (verboseLevel>1)
    G4cout << GetProcessName() << ": Track reflected" << G4endl;

  G4ThreeVector momDir = aStep.GetPostStepPoint()->GetMomentumDirection();
  if (verboseLevel>2)
    G4cout << " Old momentum direction " << momDir << G4endl;

  // Specular reflecton reverses momentum along normal
  G4double momNorm = momDir * surfNorm;
  momDir -= 2.*momNorm*surfNorm;
  
  if (verboseLevel>2)
    G4cout << " New momentum direction " << momDir << G4endl;
  
  aParticleChange.ProposeMomentumDirection(momDir);
  return &aParticleChange;
}

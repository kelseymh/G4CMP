// $Id$
//
// 20140331  Inherit from G4CMPVDriftProcess to get subtype enforcement
// 20141029  Get output hits file from configuration manager
// 20150122  Use verboseLevel instead of compiler flag for debugging
// 20150212  Remove file IO. Use sensitive detectors instead

#include "G4CMPVDriftBoundaryProcess.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4CMPMeshElectricField.hh"
#include "G4CMPSurfaceProperty.hh"
#include "G4LatticeManager.hh"
#include "G4LatticeLogical.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"
#include "G4FieldManager.hh"
#include "G4Field.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4VParticleChange.hh"
#include "G4GeometryTolerance.hh"
#include "G4SystemOfUnits.hh"

G4CMPVDriftBoundaryProcess::G4CMPVDriftBoundaryProcess(const G4String& name,
                                         const G4ParticleDefinition* carrier)
  : G4CMPVDriftProcess("G4CMP"+name+"Boundary", fChargeBoundary),
    kCarTolerance(G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()),
    theCarrier(carrier), shortName(name) {
  if (verboseLevel) G4cout << GetProcessName() << " is created " << G4endl;
}

G4CMPVDriftBoundaryProcess::~G4CMPVDriftBoundaryProcess() {}


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
  if (postStepPoint->GetStepStatus()!=fGeomBoundary ||
      aTrack.GetStepLength()<=kCarTolerance/2.) {
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }

  if (verboseLevel) G4cout << GetProcessName() << "::PostStepDoIt" << G4endl;

  // do nothing if the current step is inbound from the original volume
  G4LatticePhysical* volLattice =
    G4LatticeManager::GetLatticeManager()->GetLattice(preStepPoint->GetPhysicalVolume());
  if (volLattice != theLattice) {
    if (verboseLevel>1) {
      G4cout << GetProcessName() << ": Track inbound after reflection"
	     << G4endl;
    }

    return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);      
  }

  // Grab surface information
  G4VPhysicalVolume* thePrePV = preStepPoint->GetPhysicalVolume();
  G4VPhysicalVolume* thePostPV = postStepPoint->GetPhysicalVolume();
  G4LogicalSurface* surface = G4LogicalBorderSurface::GetSurface(thePrePV,
                                                                 thePostPV);
  G4CMPSurfaceProperty* borderSurface;

  if (surface) {
    borderSurface = static_cast <G4CMPSurfaceProperty*> (
                                                surface->GetSurfaceProperty());
  } else {
    if (verboseLevel>1) {
      G4cerr << GetProcessName() << ": No border surface defined for "
	     << thePrePV->GetName() << " to "  << thePostPV->GetName()
	     << G4endl;
    }

    return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);
  }

  absProb = borderSurface->GetAbsProb();
  absDeltaV = borderSurface->GetAbsDeltaV();
  absMinKElec = borderSurface->GetMinKElec();
  absMinKHole = borderSurface->GetMinKHole();
  electrodeV = borderSurface->GetElectrodeV();

  if (verboseLevel>1) {
    G4cout << "    Track volume: " << aTrack.GetVolume()->GetName()
	   << "\n PreStep volume: " << thePrePV->GetName()
	   << "\nPostStep volume: " << thePostPV->GetName()
	   << G4endl;
  }

  // Test #1: There is an absProb chance to be absorbed no matter what.
  if (G4UniformRand() <= absProb) {
    if (verboseLevel>1)
      G4cout << GetProcessName() << ": Track absorbed" << G4endl;

    aParticleChange.ProposeNonIonizingEnergyDeposit(GetKineticEnergy(aTrack));
    aParticleChange.ProposeTrackStatus(fStopAndKill);
    return &aParticleChange;
  }

  // Test #2: If k is larger than the threshold for this surface.
  G4LogicalVolume* LV = thePrePV->GetLogicalVolume();
  G4ThreeVector surfNorm = LV->GetSolid()->SurfaceNormal(postStepPoint->GetPosition());

  G4double absThresh = (theCarrier == G4CMPDriftElectron::Definition()) ?
                        absMinKElec : absMinKHole;

  if (GetWaveVector(aTrack).dot(surfNorm) > absThresh) {
    if (verboseLevel>1)
      G4cout << GetProcessName() << ": Track absorbed" << G4endl;

    aParticleChange.ProposeNonIonizingEnergyDeposit(GetKineticEnergy(aTrack));
    aParticleChange.ProposeTrackStatus(fStopAndKill);
    return &aParticleChange;
  }

  // Test #3: If landed on an electrode.
  G4FieldManager* fMan = LV->GetFieldManager();
  if (fMan && fMan->DoesFieldExist()) {
    const G4CMPMeshElectricField* field = 
      dynamic_cast<const G4CMPMeshElectricField*>(fMan->GetDetectorField());
    if (field) {
      G4double posVec[4] = { 4*0. };
      GetLocalPosition(aTrack, posVec);
      G4double potential = field->GetPotential(posVec);

      if((fabs(surfNorm.getZ())>0.5
            && fabs(potential-electrodeV) <= absDeltaV)) {
        if (verboseLevel>1)
          G4cout << GetProcessName() << ": Track hit electrode" << G4endl;

        aParticleChange.ProposeNonIonizingEnergyDeposit(GetKineticEnergy(aTrack));
        aParticleChange.ProposeTrackStatus(fStopAndKill);
        return &aParticleChange;
      }
    }
  } else {
    G4cout << "WTF- no field?" << G4endl;
  }

  // No absorption means reflection. Naive approach, only reflect outbound
  G4ThreeVector momentumDir = aTrack.GetMomentumDirection();
  G4double momNorm = surfNorm.dot(momentumDir);
  if (momNorm > 0.) {
    if (verboseLevel>1)
      G4cout << GetProcessName() << ": Track reflected" << G4endl;

    momentumDir -= 2.*momNorm*surfNorm;		// Simple specular reflection

    aParticleChange.ProposeMomentumDirection(momentumDir);
    aParticleChange.ProposeTrackStatus(fAlive);
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
}

G4bool G4CMPVDriftBoundaryProcess::IsApplicable(const G4ParticleDefinition& aPD) {
  return (&aPD==theCarrier);
}

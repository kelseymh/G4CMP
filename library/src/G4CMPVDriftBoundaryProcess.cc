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
#include "G4LatticeManager.hh"
#include "G4LatticeLogical.hh"
#include "G4LatticePhysical.hh"
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
  : G4CMPVDriftProcess("G4CMP"+name+"BoundaryProcess", fChargeBoundary),
    kCarTolerance(G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()),
    theCarrier(carrier), shortName(name) {
  if (verboseLevel) G4cout << GetProcessName() << " is created " << G4endl;
}

G4CMPVDriftBoundaryProcess::~G4CMPVDriftBoundaryProcess() {}


G4double 
G4CMPVDriftBoundaryProcess::GetMeanFreePath(const G4Track& /*aTrack*/,
					   G4double /*previousStepSize*/,
					   G4ForceCondition* condition) {
  *condition = Forced;
  return DBL_MAX;
}

G4VParticleChange* 
G4CMPVDriftBoundaryProcess::PostStepDoIt(const G4Track& aTrack,
					const G4Step& aStep) {    
  aParticleChange.Initialize(aTrack);
  G4StepPoint* postStepPoint = aStep.GetPostStepPoint();

  // do nothing but return if the current step is not limited by a volume
  // boundary
  if (postStepPoint->GetStepStatus()!=fGeomBoundary) { 
    return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);      
  }

  if (verboseLevel) {
    G4cout << "G4CMPDriftBoundaryProcess::PostStepDoIt" << G4endl;
  }

  // Test #1: There is an absProb chance to be absorbed no matter what.
  if (G4UniformRand() <= absProb) {
    aParticleChange.ProposeNonIonizingEnergyDeposit(GetKineticEnergy(aTrack));
    aParticleChange.ProposeTrackStatus(fStopAndKill);
    return &aParticleChange;
  }

  // Test #2: If k is larger than the threshold for this surface.
  G4VPhysicalVolume* PV = aTrack.GetVolume();
  G4ThreeVector surfNorm = PV->GetLogicalVolume()->GetSolid()->SurfaceNormal(postStepPoint->GetPosition());

  G4double absThresh;
  if (surfNorm.getZ() > 0.5)
    absThresh = absTopMinK;
  else if (surfNorm.getZ() < -0.5)
    absThresh = absBotMinK;
  else
    absThresh = absWallMinK;

  if (GetWaveVector(aTrack).dot(surfNorm) > absThresh) {
    aParticleChange.ProposeNonIonizingEnergyDeposit(GetKineticEnergy(aTrack));
    aParticleChange.ProposeTrackStatus(fStopAndKill);
    return &aParticleChange;
  }

  // Test #3: If landed on an electrode.
  G4FieldManager* fMan = PV->GetLogicalVolume()->GetFieldManager();

  if(fMan && fMan->DoesFieldExist()) {
    if(const G4CMPMeshElectricField* field = dynamic_cast<const G4CMPMeshElectricField*>(fMan->GetDetectorField())) {
      G4double posVec[4] = { 4*0. };
      GetLocalPosition(aTrack, posVec);
      G4double potential = field->GetPotential(posVec);

      if(surfNorm.getZ() > 0.5 && fabs(potential - topElectrodeV) <= absDeltaV ||
         surfNorm.getZ() < -0.5 && fabs(potential - botElectrodeV) <= absDeltaV) {
        aParticleChange.ProposeNonIonizingEnergyDeposit(GetKineticEnergy(aTrack));
        aParticleChange.ProposeTrackStatus(fStopAndKill);
        return &aParticleChange;
      }
    }
  } else {
  G4cout << "WTF- no field?" << G4endl;
  }

  // No absorption means reflection. Naive approach.
  G4ThreeVector momentumDir = aTrack.GetMomentumDirection();
  if (surfNorm.getZ() > 0.5 || surfNorm.getZ() < -0.5) {
    momentumDir.setZ(-momentumDir.getZ());
  } else {
    momentumDir.setX(-momentumDir.getX());
    momentumDir.setY(-momentumDir.getY());
  }
  aParticleChange.ProposeMomentumDirection(momentumDir);
  return &aParticleChange;
}

void G4CMPVDriftBoundaryProcess::LoadDataForTrack(const G4Track* track) {
  if (track->GetDefinition() != theCarrier) {
    G4cerr << GetProcessName() << " ERROR:  Track type "
     << track->GetDefinition()->GetParticleName() << " not valid" << G4endl;
    return;
  }

  G4CMPVDriftProcess::LoadDataForTrack(track);

  absProb = theLattice->GetAbsProb();
  absDeltaV = theLattice->GetAbsDeltaV();
  botElectrodeV = theLattice->GetBotElectrodeV();
  topElectrodeV = theLattice->GetTopElectrodeV();
  if (theCarrier == G4CMPDriftHole::Definition()) {
    absTopMinK = theLattice->GetAbsTopMinKHole();
    absBotMinK = theLattice->GetAbsBotMinKHole();
    absWallMinK = theLattice->GetAbsWallMinKHole();
  } else if (theCarrier == G4CMPDriftElectron::Definition()) {
    absTopMinK = theLattice->GetAbsTopMinKElec();
    absBotMinK = theLattice->GetAbsBotMinKElec();
    absWallMinK = theLattice->GetAbsWallMinKElec();
  }
}

G4bool G4CMPVDriftBoundaryProcess::IsApplicable(const G4ParticleDefinition& aPD) {
  return (&aPD==theCarrier);
}

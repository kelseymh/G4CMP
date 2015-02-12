// $Id$
//
// 20140331  Inherit from G4CMPVDriftProcess to get subtype enforcement
// 20141029  Get output hits file from configuration manager
// 20150122  Use verboseLevel instead of compiler flag for debugging

#include "G4CMPDriftBoundaryProcess.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4CMPConfigManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4VParticleChange.hh"
#include "G4GeometryTolerance.hh"
#include "G4SystemOfUnits.hh"

#include <fstream>


G4CMPDriftBoundaryProcess::G4CMPDriftBoundaryProcess()
  : G4CMPVDriftProcess("G4CMPDriftBoundaryProcess", fChargeBoundary),
    kCarTolerance(G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()) {
  if (verboseLevel) G4cout << GetProcessName() << " is created " << G4endl;
}

G4CMPDriftBoundaryProcess::~G4CMPDriftBoundaryProcess() { 
  file.close();
}


G4double 
G4CMPDriftBoundaryProcess::GetMeanFreePath(const G4Track& /*aTrack*/,
					   G4double /*previousStepSize*/,
					   G4ForceCondition* condition) {
  *condition = Forced;
  return DBL_MAX;
}

G4VParticleChange* 
G4CMPDriftBoundaryProcess::PostStepDoIt(const G4Track& aTrack, 
					const G4Step& aStep) {    
  aParticleChange.Initialize(aTrack);
  G4StepPoint* postStepPoint = aStep.GetPostStepPoint();

// do nothing but return if the current step is not limited by a volume boundar
  if (postStepPoint->GetStepStatus()!=fGeomBoundary) { 
    return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);      
  }

  if (verboseLevel) {
    G4cout << "G4CMPDriftBoundaryProcess::PostStepDoIt" << G4endl;
  }

  aParticleChange.ProposeNonIonizingEnergyDeposit(aTrack.GetKineticEnergy());
  aParticleChange.ProposeTrackStatus(fStopAndKill);  

  if (!file.is_open()) {
    if (verboseLevel) {
      G4cout << " opening hits file " << G4CMPConfigManager::GetHitOutput()
	     << G4endl;
    }

    file.open(G4CMPConfigManager::GetHitOutput());
  }

  file << aTrack.GetDefinition()->GetPDGCharge() << " "
       << aTrack.GetPosition().getX()/m << " "
       << aTrack.GetPosition().getY()/m << " "
       << aTrack.GetPosition().getZ()/m << " "
       << aTrack.GetGlobalTime()/ns 
       << G4endl;
  return &aParticleChange;
}

G4bool G4CMPDriftBoundaryProcess::IsApplicable(const G4ParticleDefinition& aPD)
{
  return (&aPD==G4CMPDriftElectron::Definition() ||
	  &aPD==G4CMPDriftHole::Definition());
}


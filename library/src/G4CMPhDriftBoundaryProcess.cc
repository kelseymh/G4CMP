// $Id$

#include "G4CMPhDriftBoundaryProcess.hh"
#include "G4CMPDriftHole.hh"
#include "G4LatticeManager.hh"
#include "G4LatticeLogical.hh"
#include "G4LatticePhysical.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4VParticleChange.hh"
#include "G4GeometryTolerance.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

G4CMPhDriftBoundaryProcess::G4CMPhDriftBoundaryProcess()
  : G4CMPVDriftBoundaryProcess("hDrift", G4CMPDriftHole::Definition()) {
  if (verboseLevel) G4cout << GetProcessName() << " is created " << G4endl;
}

G4CMPhDriftBoundaryProcess::~G4CMPhDriftBoundaryProcess() {}

G4ThreeVector G4CMPhDriftBoundaryProcess::GetWaveVector(const G4Track& aTrack) const {
  return aTrack.GetStep()->GetPostStepPoint()->GetMomentum() / hbarc;
}

G4double G4CMPhDriftBoundaryProcess::GetKineticEnergy(const G4Track& aTrack) const {
  return aTrack.GetStep()->GetPostStepPoint()->GetKineticEnergy();
}

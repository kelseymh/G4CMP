// $Id$

#include "G4CMPeDriftBoundaryProcess.hh"
#include "G4CMPDriftElectron.hh"
#include "G4LatticeManager.hh"
#include "G4LatticeLogical.hh"
#include "G4LatticePhysical.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4VParticleChange.hh"
#include "G4GeometryTolerance.hh"
#include "G4SystemOfUnits.hh"

G4CMPeDriftBoundaryProcess::G4CMPeDriftBoundaryProcess()
  : G4CMPVDriftBoundaryProcess("eDrift", G4CMPDriftElectron::Definition()) {;}

G4CMPeDriftBoundaryProcess::~G4CMPeDriftBoundaryProcess() {}

G4ThreeVector G4CMPeDriftBoundaryProcess::GetWaveVector(const G4Track& aTrack) const {
  G4int iv = GetValleyIndex(aTrack);
  G4ThreeVector p_local =
    GetLocalDirection(aTrack.GetStep()->GetPostStepPoint()->GetMomentum());
  return theLattice->MapPtoK_HV(iv, p_local);
}

G4double G4CMPeDriftBoundaryProcess::GetKineticEnergy(const G4Track& aTrack) const {
  G4int iv = GetValleyIndex(aTrack);
  G4ThreeVector p_local =
    GetLocalDirection(aTrack.GetStep()->GetPostStepPoint()->GetMomentum());
  return theLattice->MapPtoEkin(iv, p_local);
}

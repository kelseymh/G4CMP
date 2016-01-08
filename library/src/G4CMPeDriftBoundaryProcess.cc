// $Id$
//
// 20150601  M. Kelsey -- Follow encapsulation of boundary process actions

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

G4ThreeVector
G4CMPeDriftBoundaryProcess::GetLocalWaveVector(const G4Track& aTrack) const {
  G4int iv = GetValleyIndex(aTrack);
  return theLattice->MapV_elToK_HV(iv, GetLocalVelocityVector(aTrack));
}

// Apply kinematic absoprtion (wave-vector at surface)

G4bool G4CMPeDriftBoundaryProcess::AbsorbTrack(const G4Step& aStep) {
  return (G4CMPVDriftBoundaryProcess::AbsorbTrack(aStep) ||
    GetLocalWaveVector(*(aStep.GetTrack()))*surfNorm > absMinKElec);
}


// Apply reflection to velocity vector, not momentum

G4bool G4CMPeDriftBoundaryProcess::ReflectTrack(const G4Step& aStep) {
  G4Track* aTrack = aStep.GetTrack();
  G4ThreeVector vel = GetGlobalVelocityVector(aTrack);
  return (vel*surfNorm > 0.);		// Velocity must be outward from volume
}

G4VParticleChange* 
G4CMPeDriftBoundaryProcess::DoReflection(const G4Step& aStep) {
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
  return &aParticleChange;
}

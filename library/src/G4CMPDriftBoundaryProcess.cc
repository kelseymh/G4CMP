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
// 20160906  Follow constness of G4CMPBoundaryUtils

#include "G4CMPDriftBoundaryProcess.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4CMPGeometryUtils.hh"
#include "G4CMPSurfaceProperty.hh"
#include "G4CMPUtils.hh"
#include "G4GeometryTolerance.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4ParallelWorldProcess.hh"
#include "G4RandomDirection.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TransportationManager.hh"
#include "G4VParticleChange.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"
#include <vector>


G4CMPDriftBoundaryProcess::G4CMPDriftBoundaryProcess(const G4String& name)
  : G4CMPVDriftProcess(name, fChargeBoundary), G4CMPBoundaryUtils(this) {;}


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
  // NOTE:  G4VProcess::SetVerboseLevel is not virtual!  Can't overlaod it
  G4CMPBoundaryUtils::SetVerboseLevel(verboseLevel);

  aParticleChange.Initialize(aTrack);
  if (!IsGoodBoundary(aStep)) return &aParticleChange;

  if (verboseLevel>1) G4cout << GetProcessName() << "::PostStepDoIt" << G4endl;

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

  ApplyBoundaryAction(aTrack, aStep, aParticleChange);
  return &aParticleChange;
}


// Decide and apply different surface actions; subclasses may override

G4bool G4CMPDriftBoundaryProcess::AbsorbTrack(const G4Track& aTrack,
                                              const G4Step& aStep) const {
  if (!G4CMPBoundaryUtils::AbsorbTrack(aTrack,aStep)) return false;

  G4double absMinK = (IsElectron(&aTrack) ? GetMaterialProperty("minKElec")
		      : IsHole(&aTrack) ? GetMaterialProperty("minKHole")
		      : -1.);

  if (absMinK < 0.) {
    G4Exception("G4CMPDriftBoundaryProcess::AbsorbTrack", "Boundary003",
                EventMustBeAborted, "Invalid particle for this process.");
  }

  G4ThreeVector kvec = GetLocalWaveVector(aTrack);

  // NOTE:  K vector above is in local coords, must use local normal
  // Must use PreStepPoint volume for transform.
  G4ThreeVector surfNorm = G4CMP::GetLocalDirection(aStep.GetPreStepPoint()->GetPhysicalVolume(),
                                                    GetSurfaceNormal(aStep));

  if (verboseLevel>2) {
    G4cout << " AbsorbTrack: local k-perp " << kvec*surfNorm
	   <<" >? absMinK " << absMinK << G4endl;
  }

  return (kvec*surfNorm > absMinK);
}


// May convert recombination into phonon

void G4CMPDriftBoundaryProcess::DoAbsorption(const G4Track& aTrack,
                                             const G4Step&,
                                             G4ParticleChange&) {
  // Charge carrier gets killed and its energy goes into phonons.
  if (verboseLevel>1) {
    G4cout << GetProcessName() << "::DoAbsorption: Track absorbed" << G4endl;
  }

  G4double eKin = GetKineticEnergy(aTrack);

  //FIXME: What does the phonon distribution look like?
  G4double eDeb = theLattice->GetDebyeEnergy();
  size_t n = std::ceil(eKin / eDeb);
  aParticleChange.SetNumberOfSecondaries(n);
  while (eKin > 0.) {
    G4double E = eKin > eDeb ? eDeb : eKin;
    eKin -= eDeb;
    G4Track* sec = CreatePhonon(G4PhononPolarization::UNKNOWN,
                                G4RandomDirection(), E,
                                aTrack.GetPosition());
    aParticleChange.AddSecondary(sec);
  }

  aParticleChange.ProposeEnergy(0.);
  aParticleChange.ProposeTrackStatus(fStopButAlive);
}


void G4CMPDriftBoundaryProcess::
DoReflection(const G4Track& aTrack, const G4Step& aStep,
	     G4ParticleChange& /*aParticleChange*/) {
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
    G4VPhysicalVolume* pv = aStep.GetPreStepPoint()->GetPhysicalVolume();
    G4CMP::RotateToLocalDirection(pv, vel);
    G4ThreeVector p = theLattice->MapV_elToP(GetCurrentValley(), vel);
    G4CMP::RotateToGlobalDirection(pv, p);

    if (verboseLevel>2) {
      G4cout << " New momentum direction " << p.unit() << G4endl;

      // SANITY CHECK:  Does new momentum get back to new velocity?
      G4ThreeVector vnew = theLattice->MapPtoV_el(GetCurrentValley(),
                                                  G4CMP::GetLocalDirection(pv, p));
      G4CMP::RotateToGlobalDirection(pv, vnew);
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
}

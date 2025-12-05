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
// 20170620  Follow interface changes in G4CMPUtils, G4CMPSecondaryUtils
// 20170802  M. Kelsey -- Replace phonon production with G4CMPEnergyPartition
// 20171215  Replace boundary-point check with CheckStepBoundary()
// 20180827  M. Kelsey -- Prevent partitioner from recomputing sampling factors
// 20210328  Modify above; compute direct-phonon sampling factor here
// 20250927  AbsorbTrack() should use '&&' to require that both conditions pass
// 20251204  G4CMP-511 -- Create parallel Lambertian reflection code for charges.

#include "G4CMPDriftBoundaryProcess.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4CMPEnergyPartition.hh"
#include "G4CMPGeometryUtils.hh"
#include "G4CMPSecondaryUtils.hh"
#include "G4CMPSurfaceProperty.hh"
#include "G4CMPUtils.hh"
#include "G4GeometryTolerance.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4ParallelWorldProcess.hh"
#include "G4RandomDirection.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4StepStatus.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VParticleChange.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"
#include <vector>


G4CMPDriftBoundaryProcess::G4CMPDriftBoundaryProcess(const G4String& name)
  : G4CMPVDriftProcess(name, fChargeBoundary), G4CMPBoundaryUtils(this),
    partitioner(new G4CMPEnergyPartition) {
  partitioner->UseDownsampling(false);		// Apply preset scaling factors
}

G4CMPDriftBoundaryProcess::~G4CMPDriftBoundaryProcess() {
  delete partitioner;
}


// Process actions

G4double G4CMPDriftBoundaryProcess::
PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
				     G4double previousStepSize,
				     G4ForceCondition* condition) {
  return GetMeanFreePath(aTrack, previousStepSize, condition);
}

G4double G4CMPDriftBoundaryProcess::
GetMeanFreePath(const G4Track&, G4double, G4ForceCondition* condition) {
  *condition = Forced;
  return DBL_MAX;
}


G4VParticleChange* 
G4CMPDriftBoundaryProcess::PostStepDoIt(const G4Track& aTrack,
                                         const G4Step& aStep) {
  // NOTE:  G4VProcess::SetVerboseLevel is not virtual!  Can't overload it
  G4CMPBoundaryUtils::SetVerboseLevel(verboseLevel);

  aParticleChange.Initialize(aTrack);
  if (!IsGoodBoundary(aStep))
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);

  if (verboseLevel>1) G4cout << GetProcessName() << "::PostStepDoIt" << G4endl;

  if (verboseLevel>2) {
    if (IsElectron()) {
      G4cout << " K (" << GetValleyIndex(aTrack) << ") direction: "
	     << theLattice->MapPtoK(GetValleyIndex(aTrack),
				   GetLocalMomentum(aTrack)).unit()
	     << G4endl;
    }
    G4cout << " K direction: " << GetLocalWaveVector(aTrack).unit()
           << "\n P direction: " << GetLocalMomentum(aTrack).unit() << G4endl;
  }

  ApplyBoundaryAction(aTrack, aStep, aParticleChange);

  ClearNumberOfInteractionLengthLeft();		// All processes should do this!
  return &aParticleChange;
}


// Decide and apply different surface actions; subclasses may override

G4bool G4CMPDriftBoundaryProcess::AbsorbTrack(const G4Track& aTrack,
                                              const G4Step& aStep) const {
  G4double absMinK = (G4CMP::IsElectron(aTrack) ? GetMaterialProperty("minKElec")
		      : G4CMP::IsHole(aTrack) ? GetMaterialProperty("minKHole")
		      : -1.);

  if (absMinK < 0.) {
    G4Exception("G4CMPDriftBoundaryProcess::AbsorbTrack", "Boundary003",
                EventMustBeAborted, "Invalid particle for this process.");
  }

  G4ThreeVector kvec = GetLocalWaveVector(aTrack);

  // NOTE:  K vector above is in local coords, must use local normal
  // Must use PreStepPoint volume for transform.
  G4ThreeVector surfNorm = G4CMP::GetLocalDirection(aTrack.GetTouchable(),
                                                    G4CMP::GetSurfaceNormal(aStep));

  if (verboseLevel>2) {
    G4cout << " AbsorbTrack: local k-perp " << kvec*surfNorm
	   <<" >? absMinK " << absMinK << G4endl;
  }

  return ( (kvec*surfNorm > absMinK) &&
	   G4CMPBoundaryUtils::AbsorbTrack(aTrack, aStep) );
}


// Recombination (bandgap energy) is handled in separate AtRest process

void G4CMPDriftBoundaryProcess::DoAbsorption(const G4Track& aTrack,
                                             const G4Step&, G4ParticleChange&) {
  // Charge carrier gets killed and its energy goes into phonons.
  if (verboseLevel>1) {
    G4cout << GetProcessName() << "::DoAbsorption: Track absorbed" << G4endl;
  }

  *(G4CMPProcessUtils*)partitioner = *(G4CMPProcessUtils*)this;
  partitioner->SetVerboseLevel(verboseLevel);
  partitioner->UseVolume(aTrack.GetVolume());

  G4double eAbs = GetKineticEnergy(aTrack);

  // Compute direct-phonon downsampling here
  partitioner->ComputePhononSampling(eAbs);
  partitioner->DoPartition(0., eAbs);
  partitioner->GetSecondaries(&aParticleChange);

  if (aParticleChange.GetNumberOfSecondaries() == 0) {	// Record energy release
    aParticleChange.ProposeNonIonizingEnergyDeposit(eAbs);
  }

  aParticleChange.ProposeEnergy(0.);
  aParticleChange.ProposeTrackStatus(fStopButAlive);
}


void G4CMPDriftBoundaryProcess::
DoReflection(const G4Track& aTrack, const G4Step& aStep,
	     G4ParticleChange& /*aParticleChange*/) {
  if (verboseLevel>1)
    G4cout << GetProcessName() << ": Track reflected" << G4endl;

  // Electrons and holes need to be handled separately until we further
  // generalize the physics.

  if (IsElectron())  DoReflectionElectron(aTrack, aStep, aParticleChange);
  else if (IsHole()) DoReflectionHole(aTrack, aStep, aParticleChange);
  else {
    G4Exception("G4CMPDriftBoundaryProcess::DoReflection", "Boundary004",
                EventMustBeAborted, "Invalid particle for this process.");
  }
}

void G4CMPDriftBoundaryProcess::
DoReflectionElectron(const G4Track& aTrack, const G4Step& aStep,
		     G4ParticleChange& particleChange) {
  if (verboseLevel>1)
    G4cout << GetProcessName() << ": Electron reflected" << G4endl;

  // Get outward normal from current volume
  G4ThreeVector surfNorm = G4CMP::GetSurfaceNormal(aStep);

  // Check whether step has proper boundary-stopped geometry
  G4ThreeVector surfacePoint;
  if (!CheckStepBoundary(aStep, surfacePoint)) {
    if (verboseLevel>2)
      G4cout << " Boundary point moved to " << surfacePoint << G4endl;

    particleChange.ProposePosition(surfacePoint);	// IS THIS CORRECT?!?
  }

  // FUTURE: Get specular vs. diffuse probability from parameters
  G4bool specular = false;

  G4ThreeVector reflDir;
  if (specular) {
    G4ThreeVector vel = GetGlobalVelocityVector(aTrack);
    reflDir = DoSpecularElectron(vel, surfNorm, surfacePoint);
  } else {
    reflDir = DoDiffuseElectron(surfNorm, surfacePoint);
  }

  FillParticleChange(GetCurrentValley(), aTrack.GetKineticEnergy(), reflDir);
}

G4ThreeVector G4CMPDriftBoundaryProcess::
DoSpecularElectron(const G4ThreeVector& inDir,
		   const G4ThreeVector& surfNorm,
		   const G4ThreeVector& /*surfPos*/) const {
  if (verboseLevel>2) G4cout << " DoSpecularElectron " << surfNorm << G4endl;

  if (verboseLevel>2) {
    G4cout << " Old velocity direction " << inDir.unit() << G4endl;
  }

  // Specular reflection reverses velocity along normal
  G4double dirNorm = inDir * surfNorm;
  G4ThreeVector outDir = inDir - 2.*dirNorm*surfNorm;
  
  if (verboseLevel>2)
    G4cout << " New velocity direction " << outDir.unit() << G4endl;
  
  // Convert velocity back to momentum and update direction
  RotateToLocalDirection(outDir);
  G4ThreeVector p = theLattice->MapV_elToP(GetCurrentValley(), outDir);
  RotateToGlobalDirection(p);
  
  if (verboseLevel>2) {
    G4cout << " New momentum direction " << p.unit() << G4endl;
    
    // SANITY CHECK:  Does new momentum get back to new velocity?
    G4ThreeVector vnew = theLattice->MapPtoV_el(GetCurrentValley(),
						GetLocalDirection(p));
    RotateToGlobalDirection(vnew);
    G4cout << " Cross-check new v dir  " << vnew.unit() << G4endl;
  }

  return outDir;
}

G4ThreeVector G4CMPDriftBoundaryProcess::
DoDiffuseElectron(const G4ThreeVector& surfNorm,
		  const G4ThreeVector& /*surfPos*/) const {
  if (verboseLevel>2) G4cout << " DoDiffuseElectron " << surfNorm << G4endl;

  // Charge scatters randomly off of surface
  G4ThreeVector p = G4CMP::GetLambertianVector(surfNorm);
  return p;
}


void G4CMPDriftBoundaryProcess::
DoReflectionHole(const G4Track& /*aTrack*/, const G4Step& aStep,
		 G4ParticleChange& /*aParticleChange*/) {
  if (verboseLevel>1)
    G4cout << GetProcessName() << ": Hole reflected" << G4endl;

  G4ThreeVector surfNorm = G4CMP::GetSurfaceNormal(aStep);

  // TODO: If we do the electrons Lambertian, we should do the holes also
  G4ThreeVector momDir = aStep.GetPostStepPoint()->GetMomentumDirection();
  if (verboseLevel>2)
    G4cout << " Old momentum direction " << momDir << G4endl;
  
  G4double momNorm = momDir * surfNorm;
  momDir -= 2.*momNorm*surfNorm;
  
  if (verboseLevel>2)
    G4cout << " New momentum direction " << momDir << G4endl;
  
  aParticleChange.ProposeMomentumDirection(momDir);
}

// Called when maximum bounces have been recorded; does recombination

void G4CMPDriftBoundaryProcess::
DoFinalReflection(const G4Track& aTrack,const G4Step& aStep,
		  G4ParticleChange& aParticleChange) {
  DoAbsorption(aTrack, aStep, aParticleChange);
}

// Generate cos(theta) law for diffuse reflection, ensuring that computed
// vector is directed inward with respect to the surface normal.

G4ThreeVector
G4CMPDriftBoundaryProcess::LambertianReflection(const G4LatticePhysical* theLattice,
			   const G4ThreeVector& surfNorm, G4int valley) {
  const G4ThreeVector surfPoint = GetCurrentTrack()->GetPosition();
  return LambertianReflection(theLattice, surfNorm, valley, surfPoint);
}

G4ThreeVector
G4CMPDriftBoundaryProcess::LambertianReflection(const G4LatticePhysical* theLattice,
			   const G4ThreeVector& surfNorm, G4int valley,
			   const G4ThreeVector& surfPoint) {
  G4ThreeVector reflectedDir;
  const G4int maxTries = 1000;
  G4int nTries = 0;
  do {
    reflectedDir = G4CMP::GetLambertianVector(surfNorm);
  } while (nTries++ < maxTries &&
           !ChargeVelocityIsInward(theLattice, valley, reflectedDir, surfNorm,
                                   surfPoint));

  return reflectedDir;
}


// Check that charge is properly directed from the volume surface
// waveVector and surfNorm need to be in global coordinates

G4bool G4CMPDriftBoundaryProcess::ChargeVelocityIsInward(const G4LatticePhysical* lattice,
                                     G4int valley,
                                     const G4ThreeVector& waveVector,
                                     const G4ThreeVector& surfNorm) {
  const G4ThreeVector surfacePos = GetCurrentTrack()->GetPosition();
  return ChargeVelocityIsInward(lattice, valley, waveVector, surfNorm, surfacePos);
}

G4bool G4CMPDriftBoundaryProcess::ChargeVelocityIsInward(const G4LatticePhysical* lattice,
                                     G4int valley,
                                     const G4ThreeVector& waveVector,
                                     const G4ThreeVector& surfNorm,
                                     const G4ThreeVector& surfacePos) {
  // Get touchable for coordinate rotations
  const G4VTouchable* touchable = GetCurrentTouchable();

  if (!touchable) {
    G4Exception("ChargeVelocityIsInward", "G4CMPUtils001",
		EventMustBeAborted, "Current track does not have valid touchable!");
    return false;
  }

  // MapKtoVDir requires local direction for the wavevector
  G4ThreeVector vDir = (lattice->MapV_elToP(valley, G4CMP::GetLocalDirection(touchable, waveVector))).unit();
  // MapKToVDir is specific to phonons to vconvert wavevector to group velocity.
  // The function MapVelToP converts electron velocity vector to the momentum direction.
  // (Might need to add a .unit at the end of vDir to turn it into a unit vector)
  // The only thing here that is particle type dependent is MapKToVDir and MapVelToP.

  // Project a 1 nm step in the new direction, see if it
  // is still in the correct volume.
  G4ThreeVector localPos = G4CMP::GetLocalPosition(touchable, surfacePos);
  G4VSolid* solid = touchable->GetSolid();
  EInside trialStep = solid->Inside(localPos + 1*nm * vDir); // Is 1nm big enough? Check with Nolan.

  // Compare group velocity and surface normal in global coordinates
  G4CMP::RotateToGlobalDirection(touchable, vDir);
  return (vDir.dot(surfNorm) < 0.0 && trialStep == kInside);
}

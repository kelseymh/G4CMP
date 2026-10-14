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
// 20251015  Resolve shadowed declaration in DoFinalReflection()
// 20251024  G4CMP-519: Protect against possible zero energy in DoAbsorption()
// 20251028  G4CMP-527: Move CheckStepBoundary() to ApplyBoundaryAction()
// 20251204  G4CMP-511 -- Create parallel Lambertian reflection code for charges.
// 20251210  G4CMP-518 -- Make PhononVelocityIsInward() generic.
// 20250219  G4CMP-513 : Provide separate specular and diffuse reflection for charges.
// 20260430  W. Lamberson -- Update debugging printouts.
// 20260507  G4CMP-603 : Use specular reflection material property.

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
GetMeanFreePath(const G4Track& aTrack, G4double, G4ForceCondition* condition) {
  UpdateMeanFreePathForLatticeChangeover(aTrack);
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
  G4ThreeVector vDir = aStep.GetPreStepPoint()->GetMomentumDirection();

  // NOTE:  K vector above is in local coords, must use local normal
  // Must use PreStepPoint volume for transform.
  G4ThreeVector surfNorm = G4CMP::GetLocalDirection(aTrack.GetTouchable(),
                                                    G4CMP::GetSurfaceNormal(aStep,vDir));

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
  if (eAbs > 0.) {
    // Compute direct-phonon downsampling here
    partitioner->ComputePhononSampling(eAbs);
    partitioner->DoPartition(0., eAbs);
    partitioner->GetSecondaries(&aParticleChange);
    
    if (aParticleChange.GetNumberOfSecondaries() == 0) {
      aParticleChange.ProposeNonIonizingEnergyDeposit(eAbs);
    }
  }

  aParticleChange.ProposeEnergy(0.);
  aParticleChange.ProposeTrackStatus(fStopButAlive);
}


void G4CMPDriftBoundaryProcess::
DoReflection(const G4Track& aTrack, const G4Step& aStep,
	     G4ParticleChange& /*aParticleChange*/) {
  if (verboseLevel)
    G4cout << GetProcessName() << ": Charge reflected" << G4endl;

  G4double specProb = GetMaterialProperty("specProb");
  if (verboseLevel>2) {
    G4cout << " using " << specProb << " specular " << 1.-specProb
	   << " diffuse" << G4endl;
  }

  G4ThreeVector reflP;		// Acquire reflected momentum below

  if (G4UniformRand() < specProb) {
    reflP = DoSpecularReflection(aTrack, aStep);
  } else {			// Do diffuse reflection (electrons & holes)
    reflP = DoDiffuseReflection(aTrack, aStep);
  }

  if (verboseLevel>2)
    G4cout << " reflected momentum " << reflP << " " << reflP.mag() << G4endl;

  // Ensure that energy (not momentum!) is conserved
  FillParticleChange(GetCurrentValley(), GetKineticEnergy(aTrack), reflP);

  if (verboseLevel>3) {
    aParticleChange.DumpInfo();
    G4ThreeVector pcp =
      aParticleChange.CalcMomentum(aParticleChange.GetEnergy(),
				   *aParticleChange.GetMomentumDirection(),
				   aParticleChange.GetMass());

    G4cout << " particleChange momentum " << pcp << G4endl;
  }
}

// NOTE:  These functions return a full momentum vector, not a unit vector

G4ThreeVector G4CMPDriftBoundaryProcess::
DoDiffuseReflection(const G4Track& aTrack, const G4Step& aStep) {
  if (verboseLevel > 1) G4cout << " DoDiffuse" << G4endl;

  if (verboseLevel>3)
    G4cout << " surfNorm " << G4CMP::GetSurfaceNormal(aStep) << G4endl;
  
  G4ThreeVector reflP =
    G4CMP::LambertianReflection(theLattice, G4CMP::GetSurfaceNormal(aStep),
				GetCurrentValley());
  if (verboseLevel>2)
    G4cout << " reflP " << reflP << " " << reflP.mag() << G4endl;

  // Rescale direction to match incident momentum
  reflP *= GetGlobalMomentum(aTrack).mag();
  
  return reflP;
}

G4ThreeVector G4CMPDriftBoundaryProcess::
DoSpecularReflection(const G4Track& aTrack, const G4Step& aStep) {
  if (verboseLevel > 1) G4cout << " DoSpecular" << G4endl;

  G4ThreeVector reflP;
  
  if (IsElectron())  reflP = DoSpecularElectron(aTrack, aStep);
  else if (IsHole()) reflP = DoSpecularHole(aTrack, aStep);
  else {
    G4ExceptionDescription ed;
  ed << "Unexpected particle type reaching charge drift boundary reflection: "
     << aTrack.GetParticleDefinition()->GetParticleName()
     << "\n(only electrons and holes are supported)";
  G4Exception("G4CMPDriftBoundaryProcess::DoReflection", "Boundary004",
              FatalException, ed);
  }

  return reflP;
}

G4ThreeVector G4CMPDriftBoundaryProcess::
DoSpecularElectron(const G4Track& aTrack, const G4Step& aStep) {
  if (verboseLevel > 1) G4cout << " DoSpecularElectron" << G4endl;

  G4ThreeVector kDir = GetLocalWaveVector(aTrack).unit();
  RotateToGlobalDirection(kDir);
  G4ThreeVector surfNorm = G4CMP::GetSurfaceNormal(aStep);

  if (verboseLevel > 2) {
    G4cout << " surfNorm " << surfNorm << G4endl
	   << " Old wavevector direction " << kDir << G4endl;
  }

  // Specular reflection reverses wavevector along normal
  G4double dirNorm = kDir * surfNorm;
  kDir -= 2.*dirNorm*surfNorm;

  if (verboseLevel > 2)
    G4cout << " Trying reflected wavevector " << kDir << G4endl;

  // If reflected velocity is outward facing, fall back to diffuse reflection
  // NOTE: VelocityIsInward() expects to receive _global_ wavevector
  G4int ivalley = GetCurrentValley();
  if (!G4CMP::VelocityIsInward(theLattice, ivalley, GetGlobalDirection(kDir),
			       surfNorm)) {
    if (verboseLevel>2) G4cout << " specular reflection failed" << G4endl;
    return DoDiffuseReflection(aTrack, aStep);
  }
  
  if (verboseLevel > 2) {
    G4cout << " New wavevector direction " << kDir << G4endl;
  }
  
  // Convert wavevector back to momentum and update direction
  G4ThreeVector pDir = theLattice->MapV_elToP(GetCurrentValley(),
					      GetLocalDirection(kDir)).unit();
  RotateToGlobalDirection(pDir);
  pDir *= GetLocalMomentum(aTrack).mag();

  return pDir;
}

G4ThreeVector G4CMPDriftBoundaryProcess::
DoSpecularHole(const G4Track& aTrack, const G4Step& aStep) {
  if (verboseLevel>1) G4cout << " DoSpecularHole" << G4endl;

  G4ThreeVector surfNorm = G4CMP::GetSurfaceNormal(aStep);

  G4ThreeVector momDir = aStep.GetPostStepPoint()->GetMomentumDirection();
  if (verboseLevel > 2) {
    G4cout << " Old momentum direction " << momDir << G4endl;
  }
  
  G4double momNorm = momDir * surfNorm;
  momDir -= 2.*momNorm*surfNorm;
  
  if (verboseLevel>2) {
    G4cout << " New momentum direction " << momDir << G4endl;
  }
  
  return aTrack.GetMomentum().mag() * momDir;
}


// Called when maximum bounces have been recorded; does recombination

void G4CMPDriftBoundaryProcess::
DoFinalReflection(const G4Track& aTrack,const G4Step& aStep,
		  G4ParticleChange& particleChange) {
  DoAbsorption(aTrack, aStep, particleChange);
}

/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPBoundaryUtils.cc
/// \brief Implementation of the G4CMPBoundaryUtils class
///   Provides useful general functions for use by boundary processes,
///   which do not have a common base class for inheritance.
///
///   Use via multiple inheritance with concrete boundary classes.
//
// $Id$
//
// 20160904  Add electrode pattern handling
// 20160906  Make most functions const, provide casting function for matTable
// 20161114  Use G4CMPVTrackInfo
// 20170710  Look for skin surface (LV) if border surface not found
// 20170713  Report undefined surfaces only once per job, not a failure
// 20171215  Change 'CheckStepStatus()' to 'IsBoundaryStep()', add function
//	     to validate step trajectory to boundary.
// 20201112  Add warning message to base DoTransmission() function (c.f.
//	     warning message in base DoReflection()).  Pass verbosity through
//	     to electrode.
// 20210923  Use >= in maximum reflections check.
// 20211207  Replace G4Logical*Surface with G4CMP-specific versions.
// 20250413  Protect debugging messages with verbosity.
// 20250415  Suppress error for same PV if starting at boundary.
// 20250423  Remove error suppression for starting at boundary.

#include "G4CMPBoundaryUtils.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPGeometryUtils.hh"
#include "G4CMPLogicalBorderSurface.hh"
#include "G4CMPLogicalSkinSurface.hh"
#include "G4CMPSurfaceProperty.hh"
#include "G4CMPProcessUtils.hh"
#include "G4CMPVTrackInfo.hh"
#include "G4CMPTrackUtils.hh"
#include "G4CMPUtils.hh"
#include "G4CMPVElectrodePattern.hh"
#include "G4ExceptionSeverity.hh"
#include "G4GeometryTolerance.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalSurface.hh"
#include "G4ParticleChange.hh"
#include "G4ParticleDefinition.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"
#include "G4VSolid.hh"


// Constructor and destructor

G4CMPBoundaryUtils::G4CMPBoundaryUtils(G4VProcess* process)
  : buVerboseLevel(G4CMPConfigManager::GetVerboseLevel()),
    procName(process->GetProcessName()), procUtils(0),
    kCarTolerance(G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()),
    maximumReflections(-1), prePV(0), postPV(0), surfProp(0), matTable(0),
    electrode(0) {
  procUtils = dynamic_cast<G4CMPProcessUtils*>(process);
  if (!procUtils) {
    G4Exception("G4CMPBoundaryUtils::G4CMPBoundaryUtils", "Boundary000",
		FatalException, "Must be passed a G4CMP process!");
  }
}

G4CMPBoundaryUtils::~G4CMPBoundaryUtils() {;}

// Initialize volumes, surface properties, etc.

G4bool G4CMPBoundaryUtils::IsGoodBoundary(const G4Step& aStep) {
  const G4ParticleDefinition* pd = aStep.GetTrack()->GetParticleDefinition();
  maximumReflections = 
    (G4CMP::IsChargeCarrier(pd) ? G4CMPConfigManager::GetMaxChargeBounces()
     : G4CMP::IsPhonon(pd) ? G4CMPConfigManager::GetMaxPhononBounces() : -1);

  if (buVerboseLevel>1) {
    G4cout << procName << "::IsGoodBoundary maxRefl " << maximumReflections
	   << G4endl;
  }

  return (IsBounaryStep(aStep) &&
	  GetBoundingVolumes(aStep) &&
	  GetSurfaceProperty(aStep));
}

G4bool G4CMPBoundaryUtils::IsBounaryStep(const G4Step& aStep) {
  if (buVerboseLevel>1) {
    G4cout << procName << "::IsBounaryStep status "
	   << aStep.GetPostStepPoint()->GetStepStatus()
	   << " length " << aStep.GetStepLength() << G4endl;
  }

  // do nothing if the current step is not limited by a volume boundary,
  // or if it is the returning "null step" after a reflection
  return aStep.GetPostStepPoint()->GetStepStatus() == fGeomBoundary;
}

G4bool G4CMPBoundaryUtils::GetBoundingVolumes(const G4Step& aStep) {
  prePV = aStep.GetPreStepPoint()->GetPhysicalVolume();
  postPV = aStep.GetPostStepPoint()->GetPhysicalVolume();

  if (prePV == postPV) {
    if (buVerboseLevel) {
      G4cerr << procName << " ERROR: fGeomBoundary status set, but"
	     << " pre- and post-step volumes are identical!" << G4endl;
    }
    return false;
  }

  // do nothing if the current step is inbound from outside the original volume
  G4LatticePhysical* volLattice =
    G4LatticeManager::GetLatticeManager()->GetLattice(prePV);
  if (volLattice != procUtils->GetLattice()) {
    if (buVerboseLevel>1) {
      G4cout << procName << ": Track inbound after reflection" << G4endl;
    }
    return false;
  }

  if (buVerboseLevel>1) {
    G4cout <<   "  PreStep volume: " << prePV->GetName() << " @ "
	   << aStep.GetPreStepPoint()->GetPosition()
	   << "\n PostStep volume: " << (postPV?postPV->GetName():"OutOfWorld")
	   << " @ " << aStep.GetPostStepPoint()->GetPosition()
	   << G4endl;
  }

  return true;
}

G4bool G4CMPBoundaryUtils::GetSurfaceProperty(const G4Step& aStep) {
  surfProp = nullptr;				// Avoid stale cache!
  matTable = nullptr;
  electrode = nullptr;
  
  // Look for specific surface between pre- and post-step points first
  G4LogicalSurface* surface =
    G4CMPLogicalBorderSurface::GetSurface(prePV, postPV);
  if (!surface) {			// Then for generic pre-setp surface
    surface = G4CMPLogicalSkinSurface::GetSurface(prePV->GetLogicalVolume());
  }

  BoundaryPV bound(prePV,postPV);	// Avoid multiple temporaries below

  // Report missing surface once per boundary
  if ((hasSurface.find(bound) == hasSurface.end())
      && !surface) {
    G4Exception((procName+"::GetSurfaceProperty").c_str(), "Boundary001",
                JustWarning, ("No surface defined between " +
			      prePV->GetName() + " and " +
			      postPV->GetName()).c_str());
  }

  hasSurface[bound] = false;		// Remember this boundary

  if (!surface) return true;			// Can handle undefined surfaces

  G4SurfaceProperty* baseSP = surface->GetSurfaceProperty();
  if (!baseSP) {
    G4Exception((procName+"::GetSurfaceProperty").c_str(),
		"Boundary002", JustWarning,
		("No surface property defined for "+surface->GetName()).c_str()
		);
    return true;			// Can handle undefined surfaces
  }

  // Verify that surface property is G4CMP compatible
  surfProp = dynamic_cast<G4CMPSurfaceProperty*>(baseSP);
  if (!surfProp) {
    G4Exception((procName+"::GetSurfaceProperty").c_str(),
		"Boundary003", EventMustBeAborted,
		"Surface property is not G4CMP compatible");
    return false;			// Badly defined, not undefined!
  }
    
  // Extract particle-specific information for later
  const G4ParticleDefinition* pd = aStep.GetTrack()->GetParticleDefinition();
  if (G4CMP::IsChargeCarrier(pd)) {
    matTable = surfProp->GetChargeMaterialPropertiesTablePointer();
    electrode = surfProp->GetChargeElectrode();
  }

  if (G4CMP::IsPhonon(pd)) {
    matTable = surfProp->GetPhononMaterialPropertiesTablePointer();
    electrode = surfProp->GetPhononElectrode();
  }

  if (!matTable) {
    G4Exception((procName+"::GetSurfaceProperty").c_str(),
		"Boundary004", JustWarning,
		(pd->GetParticleName()+" has no surface properties").c_str()
		);
    return true;			// Can handle undefined surfaces
  }

  // Initialize electrode for current track
  if (electrode) {
    electrode->SetVerboseLevel(buVerboseLevel);
    electrode->LoadDataForTrack(aStep.GetTrack());
  }

  hasSurface[bound] = true;		// Record good surface defined

  return true;
}


// Check whether end of step is actually on surface of volume
// "surfacePoint" returns post-step position, or computed surface point

G4bool G4CMPBoundaryUtils::CheckStepBoundary(const G4Step& aStep,
					     G4ThreeVector& surfacePoint) {
  G4StepPoint* preP = aStep.GetPreStepPoint();
  G4StepPoint* postP = aStep.GetPostStepPoint();
  GetBoundingVolumes(aStep);
  surfacePoint = postP->GetPosition();		// Correct if valid boundary

  // Get pre- and post-step positions in pre-step volume coordinates
  G4VSolid* preSolid = prePV->GetLogicalVolume()->GetSolid();

  G4ThreeVector prePos = preP->GetPosition();
  G4CMP::RotateToLocalPosition(preP->GetTouchable(), prePos);

  G4ThreeVector postPos = surfacePoint;
  G4CMP::RotateToLocalPosition(preP->GetTouchable(), postPos);

  if (buVerboseLevel>2) {
    G4cout << "CheckStepBoundary: in prePV (" << prePV->GetName() << ") frame"
	   << "\n  preStep @ " << prePos << "\n postStep @ " << postPos
	   << G4endl;
  }

  // Verify that post-step position is on surface of volume
  EInside postIn = preSolid->Inside(postPos);

  if (buVerboseLevel>2) {
    G4cout << "\n Is postStep location on surface of preStep Volume? "
	   << (postIn==kOutside ? "outside" :
	       postIn==kInside  ? "inside" :
	       postIn==kSurface ? "surface" : "INVALID") << G4endl;
  }

  // If post-step position not proper surface point, compute intersection
  if (postIn != kSurface) {
    if (buVerboseLevel>2)
      G4cout << " OLD SURFACE POINT: " << surfacePoint << G4endl;

    G4ThreeVector along = (postPos-prePos).unit();	// Trajectory direction
    surfacePoint = prePos + preSolid->DistanceToOut(prePos,along)*along;
    if (buVerboseLevel>2) {
      G4cout << " moving preStep by " << preSolid->DistanceToOut(prePos,along)
	     << " along " << along << G4endl
	     << " NEW SURFACE POINT: " << surfacePoint << G4endl;
    }

    postIn = preSolid->Inside(surfacePoint);
    if (buVerboseLevel>2) {
      G4cout << "\n Is adjusted location on surface of preStep Volume? "
	     << (postIn==kOutside ? "outside" :
		 postIn==kInside  ? "inside" :
		 postIn==kSurface ? "surface" : "INVALID") << G4endl;
    }

    // Double check calculation -- point "must" now be on surface!
    if (postIn != kSurface) {
      G4Exception((procName+"::CheckBoundaryPoint").c_str(),
		  "Boundary005", EventMustBeAborted,
		  "Boundary-limited step cannot find boundary surface point"
		  );
      return false;
    }

    // Move surface point to world coordinate system
    G4CMP::RotateToGlobalPosition(preP->GetTouchable(), surfacePoint);
  }

  return (postIn == kSurface);
}


// Implement PostStepDoIt() in a common way; processes should call through

void 
G4CMPBoundaryUtils::ApplyBoundaryAction(const G4Track& aTrack,
					const G4Step& aStep,
					G4ParticleChange& aParticleChange) {
  aParticleChange.Initialize(aTrack);

  if (!matTable) {
    if (buVerboseLevel>2) G4cout << "BU::Apply: !matTable" << G4endl;
    DoSimpleKill(aTrack, aStep, aParticleChange);
  } else if (electrode && electrode->IsNearElectrode(aStep)) {
    if (buVerboseLevel>2) G4cout << "BU::Apply: absorb at electrode" << G4endl;
    electrode->AbsorbAtElectrode(aTrack, aStep, aParticleChange);
  } else if (AbsorbTrack(aTrack, aStep)) {
    if (buVerboseLevel>2) G4cout << "BU::Apply: Absorption" << G4endl;
    DoAbsorption(aTrack, aStep, aParticleChange);
  } else if (MaximumReflections(aTrack)) {
    if (buVerboseLevel>2) G4cout << "BU::Apply: maxRef" << G4endl;
    DoSimpleKill(aTrack, aStep, aParticleChange);
  } else if (ReflectTrack(aTrack, aStep)) {
    if (buVerboseLevel>2) G4cout << "BU::Apply: Reflection" << G4endl;
    DoReflection(aTrack, aStep, aParticleChange);
  } else {
    if (buVerboseLevel>2) G4cout << "BU::Apply: Transmission" << G4endl;
    DoTransmission(aTrack, aStep, aParticleChange);
  }
}


// Default conditions for absorption or reflection

G4bool G4CMPBoundaryUtils::AbsorbTrack(const G4Track&, const G4Step&) const {
  G4double absProb = GetMaterialProperty("absProb");
  if (buVerboseLevel>2)
    G4cout << " AbsorbTrack: absProb " << absProb << G4endl;

  return (G4UniformRand() <= absProb);
}

G4bool G4CMPBoundaryUtils::ReflectTrack(const G4Track&, const G4Step&) const {
  G4double reflProb = GetMaterialProperty("reflProb");
  if (buVerboseLevel>2)
    G4cout << " ReflectTrack: reflProb " << reflProb << G4endl;

  return (G4UniformRand() <= reflProb);
}

G4bool G4CMPBoundaryUtils::MaximumReflections(const G4Track& aTrack) const {
  auto trackInfo = G4CMP::GetTrackInfo<G4CMPVTrackInfo>(aTrack);
  trackInfo->IncrementReflectionCount();

  return (maximumReflections >= 0 &&
    trackInfo->ReflectionCount() >= static_cast<size_t>(maximumReflections));
}


// Simple absorption deposits non-ionizing energy

void G4CMPBoundaryUtils::DoAbsorption(const G4Track& aTrack,
				      const G4Step& /*aStep*/,
				      G4ParticleChange& aParticleChange) {
  if (buVerboseLevel>1) G4cout << procName << ": Track absorbed" << G4endl;

  G4double ekin = procUtils->GetKineticEnergy(aTrack);
  aParticleChange.ProposeNonIonizingEnergyDeposit(ekin);
  aParticleChange.ProposeTrackStatus(fStopAndKill);
  aParticleChange.ProposeEnergy(0.);
}

void G4CMPBoundaryUtils::DoReflection(const G4Track& aTrack,
				      const G4Step& aStep,
				      G4ParticleChange& aParticleChange) {
  G4cerr << procName << " WARNING!  G4CMPBoundaryUtils::DoReflection invoked."
	 << "\n Process should have overridden this version!"
	 << "  Results may be non-physical" << G4endl;

  if (buVerboseLevel>1) {
    G4cout << procName << ": Track reflected "
           << G4CMP::GetTrackInfo<G4CMPVTrackInfo>(aTrack)->ReflectionCount()
	   << " times." << G4endl;
  }

  G4ThreeVector pdir = aTrack.GetMomentumDirection();
  G4ThreeVector norm = G4CMP::GetSurfaceNormal(aStep);	// Outward normal
  pdir -= 2.*(pdir.dot(norm))*norm;			// Reverse along normal

  aParticleChange.ProposeMomentumDirection(pdir);
}

void G4CMPBoundaryUtils::DoSimpleKill(const G4Track& /*aTrack*/,
				      const G4Step& /*aStep*/,
				      G4ParticleChange& aParticleChange) {
  if (buVerboseLevel>1) G4cout << procName << ": Track killed" << G4endl;

  aParticleChange.ProposeTrackStatus(fStopAndKill);
}

void 
G4CMPBoundaryUtils::DoTransmission(const G4Track& aTrack,
				   const G4Step& aStep,
				   G4ParticleChange& aParticleChange) {
  G4cerr << procName << " WARNING!  G4CMPBoundaryUtils::DoTransmission invoked."
	 << "\n Process should have overridden this version!"
	 << "  Track will be killed as leaving volume" << G4endl;

  if (buVerboseLevel>1)
    G4cout << procName << ": Track transmission requested" << G4endl;

  DoSimpleKill(aTrack, aStep, aParticleChange);
}


// Access information from materials table even when const

G4double G4CMPBoundaryUtils::GetMaterialProperty(const G4String& key) const {
  return const_cast<G4MaterialPropertiesTable*>(matTable)->GetConstProperty(key);
}

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

#include "G4CMPBoundaryUtils.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPSurfaceProperty.hh"
#include "G4CMPProcessUtils.hh"
#include "G4CMPTrackInformation.hh"
#include "G4CMPUtils.hh"
#include "G4CMPVElectrodePattern.hh"
#include "G4ExceptionSeverity.hh"
#include "G4GeometryTolerance.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSurface.hh"
#include "G4Navigator.hh"
#include "G4ParticleChange.hh"
#include "G4ParticleDefinition.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"


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
    G4cout << procName << "::LoadDataForStep maxRefl " << maximumReflections
	   << G4endl;
  }

  return (CheckStepStatus(aStep) &&
	  GetBoundingVolumes(aStep) &&
	  GetSurfaceProperty(aStep));
}

G4bool G4CMPBoundaryUtils::CheckStepStatus(const G4Step& aStep) {
  if (buVerboseLevel>1) {
    G4cout << procName << "::CheckStepStatus status "
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
	   << "\n PostStep volume: " << postPV->GetName() << " @ "
	   << aStep.GetPostStepPoint()->GetPosition()
	   << G4endl;
  }

  return true;
}

G4bool G4CMPBoundaryUtils::GetSurfaceProperty(const G4Step& aStep) {
  G4LogicalSurface* surface = G4LogicalBorderSurface::GetSurface(prePV, postPV);
  if (!surface) {
    G4Exception((procName+"::GetSurfaceProperty").c_str(), "Boundary001",
                EventMustBeAborted, ("No surface defined between " +
                                    prePV->GetName() + " and " +
                                    postPV->GetName() + ".").c_str());
    return false;
  }

  G4SurfaceProperty* baseSP = surface->GetSurfaceProperty();
  if (!baseSP) {
    G4Exception((procName+"::GetSurfaceProperty").c_str(),
		"Boundary002", EventMustBeAborted,
		("No surface property defined for "+surface->GetName()).c_str()
		);
    return false;
  }

  // Verify that surface property is G4CMP compatible
  surfProp = dynamic_cast<G4CMPSurfaceProperty*>(baseSP);
  if (!surfProp) {
    G4Exception((procName+"::GetSurfaceProperty").c_str(),
		"Boundary003", EventMustBeAborted,
		"Surface property is not G4CMP compatible");
    return false;
  }
    
  // Extract particle-specific information for later
  const G4ParticleDefinition* pd = aStep.GetTrack()->GetParticleDefinition();
  const G4MaterialPropertiesTable* constTbl = 0;
  if (G4CMP::IsChargeCarrier(pd)) {
    constTbl = surfProp->GetChargeMaterialPropertiesTablePointer();
    electrode = surfProp->GetChargeElectrode();
  }

  if (G4CMP::IsPhonon(pd)) {
    constTbl = surfProp->GetPhononMaterialPropertiesTablePointer();
    electrode = surfProp->GetPhononElectrode();
  }

  if (!constTbl) {
    G4Exception((procName+"::GetSurfaceProperty").c_str(),
		"Boundary004", EventMustBeAborted,
		(pd->GetParticleName()+" has no material properties").c_str()
		);
    return false;
  }

  // Must store non-const pointer, as table has no const accessors
  matTable = const_cast<G4MaterialPropertiesTable*>(constTbl);

  return true;
}


// Implement PostStepDoIt() in a common way; processes should call through

void G4CMPBoundaryUtils::
ApplyBoundaryAction(const G4Track& aTrack, const G4Step& aStep,
		    G4ParticleChange& aParticleChange) {
  aParticleChange.Initialize(aTrack);

  if (!IsGoodBoundary(aStep)) return;		// May have been done already

  // If the particle doesn't get absorbed, it either reflects or transmits
  if (electrode && electrode->IsNearElectrode(aStep)) {
    if (buVerboseLevel>1) G4cout << procName << ": Track electrode " << G4endl;
    electrode->AbsorbAtElectrode(aTrack, aStep, aParticleChange);
  } else if (AbsorbTrack(aTrack, aStep)) {
    if (electrode) {
      DoSimpleKill(aTrack, aStep, aParticleChange);
    } else {
      if (buVerboseLevel>1) G4cout << procName << ": Track absorbed " << G4endl;
      DoAbsorption(aTrack, aStep, aParticleChange);
    }
  } else if (MaximumReflections(aTrack)) {
    DoSimpleKill(aTrack, aStep, aParticleChange);
  } else if (ReflectTrack(aTrack, aStep)) {
    DoReflection(aTrack, aStep, aParticleChange);
  } else {
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
  G4CMPTrackInformation* trackInfo = procUtils->GetTrackInfo(aTrack);
  trackInfo->IncrementReflectionCount();

  return (maximumReflections >= 0 &&
	  trackInfo->GetReflectionCount() > maximumReflections);
}


// Simple absorption deposits non-ionizing energy

void G4CMPBoundaryUtils::DoAbsorption(const G4Track& aTrack,
				      const G4Step& /*aStep*/,
				      G4ParticleChange& aParticleChange) {
  if (buVerboseLevel>1) G4cout << procName << ": Track absorbed" << G4endl;

  G4double ekin = procUtils->GetKineticEnergy(aTrack);
  aParticleChange.ProposeNonIonizingEnergyDeposit(ekin);
  aParticleChange.ProposeTrackStatus(fStopButAlive);
}

void G4CMPBoundaryUtils::DoReflection(const G4Track& aTrack,
				      const G4Step& aStep,
				      G4ParticleChange& aParticleChange) {
  G4cerr << procName << " WARNING!  G4CMPBoundaryUtils::DoReflection invoked."
	 << "\n Process should have overridden this version!"
	 << "  Results may be non-physical" << G4endl;

  if (buVerboseLevel>1) {
    G4cout << procName << ": Track reflected "
           << procUtils->GetTrackInfo(aTrack)->GetReflectionCount()
	   << " times." << G4endl;
  }

  G4ThreeVector pdir = aTrack.GetMomentumDirection();
  G4ThreeVector norm = procUtils->GetSurfaceNormal(aStep);	// Outward normal
  pdir -= 2.*(pdir.dot(norm))*norm;			// Reverse along normal

  aParticleChange.ProposeMomentumDirection(pdir);
}

void G4CMPBoundaryUtils::DoSimpleKill(const G4Track& aTrack,
				      const G4Step& aStep,
				      G4ParticleChange& aParticleChange) {
  if (buVerboseLevel>1) G4cout << procName << ": Track killed" << G4endl;

  aParticleChange.ProposeTrackStatus(fStopButAlive);
}

void 
G4CMPBoundaryUtils::DoTransmission(const G4Track& aTrack,
				   const G4Step& aStep,
				   G4ParticleChange& aParticleChange) {
  if (buVerboseLevel>1) G4cout << procName << ": Track transmitted" << G4endl;

  DoSimpleKill(aTrack, aStep, aParticleChange);
}


// Access information from materials table even when const

G4double G4CMPBoundaryUtils::GetMaterialProperty(const G4String& key) const {
  return const_cast<G4MaterialPropertiesTable*>(matTable)->GetConstProperty(key);
}

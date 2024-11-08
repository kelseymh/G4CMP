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
     : G4CMP::IsPhonon(pd) ? G4CMPConfigManager::GetMaxPhononBounces()
     : G4CMP::IsBogoliubovQP(pd) ? G4CMPConfigManager::GetMaxBogoliubovQPBounces() : -1);

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



  //Note that since now the GetMFP function has called a lattice update, procUtils->GetLattice() should just be the "current" lattice we're in.
  //  G4cout << "REL inside the GetBoundingVolumes function. aStep.GetStepLength() = " << aStep.GetStepLength() << G4endl;
  G4cout << "---- In GetBoundingVolumes: procUtils->GetLattice(): " << procUtils->GetLattice() << ", prePVlattice: " << G4LatticeManager::GetLatticeManager()->GetLattice(prePV) << ", volLattice (postPVLattice): " << G4LatticeManager::GetLatticeManager()->GetLattice(postPV) << G4endl;
  G4LatticePhysical* volLattice = G4LatticeManager::GetLatticeManager()->GetLattice(postPV);

  //First scenario: we DON'T have a lattice in the pre-step point. This happens at least when we're on the boundary of a lattice and                          
  //the world. Here the phonon "turns around" in the world volume, in a step that is of negligible length. This negligible length is                          
  //typically around 1E-15, so to account for these safely we use a tolerance of 1E-12, which is well below the physics scales                                
  //relevant in these kinds of sims.
  double stepLengthTolerance = 1E-12 * CLHEP::m;
  if (G4LatticeManager::GetLatticeManager()->GetLattice(prePV) == 0 ){
    G4cout << "---- In GetBoundingVolumes: Current lattice is zero." << G4endl;
    
    //First: if the current (i.e. pre-step, procUtils->GetLattice()) lattice is the same as post-step volume lattice.
    //This occurs, if, for example, the current volume is World/vacuum, and a lattice changeover/update failed in the MFP step because
    //there is no lattice to update to. Hence, the "current" procUtils->GetLattice() lattice is really the prePV one from the *previous* step.
    if( volLattice == procUtils->GetLattice() ){
      G4cout << "----- In GetBoundingVolumes: Current Lattice is equal to post-PV lattice." << G4endl;
      
      //If the step length is below tolerance, we need to return false so we don't try to "double-count" the boundary action.
      //The small step sizes occur when the phonon "turns around" on a boundary with a volume that doesn't have a lattice.
      if( aStep.GetStepLength() <= stepLengthTolerance ){
        G4cout << "------ In GetBoundingVolumes: Step length is below step length tolerance." << G4endl;
        return false;
      }
      //If the step length is long. I think this can only happen if we're somehow having a long track in a non-lattice volume. For now
      //we'll keep the control block here just in case (and make it return what would be returned in its absence anyway).
      else{
        G4cout << "------ In GetBoundingVolumes: Step length is above step length tolerance." << G4endl;
        return true; //TBD, but set to be consistent with older version
      }
    }
    //If the current lattice is not equal to the post-step lattice. From a first glance at phonon dynamics I don't think this happens,
    //but we'll keep the control block here so that if we see it does happen, we can make the call then.
    else{
      G4cout << "----- In GetBoundingVolumes: Current Lattice is not equal to post-PV lattice." << G4endl;
      return true; //TBD, but set to be consistent with older version
    }
  }
  
  //Second scenario: we DO have a lattice in the pre-step point. In this case, we have two options. The initial step is one in which
  //a phonon from lattice 1 is approaching lattice 2, and has a finite step length. Here the "procUtils, current" lattice is the
  //one from which the phonon is incident, and the "postPV, volLattice" is the one that it is approaching. They are different and
  //as a result, the internal (second) if statement fires here, and as long as this first step is finite in extent, you run the
  //logic to see if you do a transmission or reflection, etc.
  // - For transmission, that's all she wrote: there is no turnaround step, and the phonon is transferred to the next lattice (see
  //   doTransmission in the phononBoundaryProcess class).  
  // - For reflection, there is an additional "infinitesimal" step, in which lattice 2 becomes the "current" lattice and lattice 1
  //   becomes the "far" lattice. The same (second) if statement triggers, but now the step length is tiny, which we can flag. For
  //   this step, we don't want to trigger any of the doTransmission/doReflection logic (as it would be "double-counting" that physics),
  //   so we return false. 
  if (G4LatticeManager::GetLatticeManager()->GetLattice(prePV) != 0 ){
    G4cout << "---- In GetBoundingVolumes: Current lattice is not zero." << G4endl;
    
    //If the current (i.e. pre-step, procUtils->GetLattice()) lattice is different from the post-step volume lattice, then this is where
    //some important logic must happen.
    if( volLattice != procUtils->GetLattice() ){
      G4cout << "----- In GetBoundingVolumes: Current Lattice ("<< procUtils->GetLattice() << ") is not equal to post-PV lattice (" << volLattice << ")." << G4endl;

      //If the step length is tiny, this means that we're in a "turnaround" step characteristic of reflection back into the lattice
      //that we approached from. Here we return false so that we don't try to run another boundary process action for this turnaround
      //step -- that has already been done
      if( aStep.GetStepLength() <= stepLengthTolerance ){
        G4cout << "------ In GetBoundingVolumes: Step length is below step length tolerance." << G4endl;
        return false;
      }
      //Otherwise, the track is actually moving through a volume before it hits this surface, and it needs to run the logic to see if
      //it reflects, transmits, etc. Need to return true so that logic can run.
      else{
        G4cout << "------ In GetBoundingVolumes: Step length is above step length tolerance." << G4endl; 
        return true;
      }
    }
    //If the current lattice is NOT different from the post-step lattice, this should have already been covered at the beginning of this
    //function -- throw an error here.
    else{
      G4cout << "\n\n\nREL Current lattice is not different from postPV lattice. This should have been covered already -- figure out why it hasn't.\n\n\n" << G4endl;
      return false; //TBD
    }
  }



  /*  
  //REL UPDATED 11/5/2024 -- this original code does not work, because when we inevitably do an update of the lattice in the MFP calculation at any given step,
  //                         we end up having the volLattice and the procUtils->GetLattice() values be the same even for the turnaround step.
  //REL 11/5/2024 -- this is the original code. It should still be applicable. The reason is twofold:
  //1. Before the phonon is transmitted across a boundary, its procUtils->GetLattice() pulls a lattice set by *track* information. As long as
  //   the phonon is in the original volume (i.e. hasn't been transmitted), the lattice here should remain the same. This holds even for when
  //   a phonon is hitting a boundary with another lattice that IS well defined. Here, the procUtils->GetLattice() should still return the
  //   original volume since we haven't used a SetLattice() function that uses the volume-based info.
  //2. During the transmission process, the procUtils->SetLattice() function is called and the lattice is updated with the new volume. Since
  //   GetBoundaryVolumes() is called twice per step but both times are BEFORE this setLattice instance, a step with a phonon incident on a
  //   boundary will not end up thinking the volLattice and procUtils->GetLattice() are different, and won't block any action. In a scenario
  //   where there's a "reflection" step that is infinitesimal, you've technically hopped over to a new lattice (so the volLattice has changed),
  //   but since there's no explicit change to the procUtils lattice, it's still the *original* lattice. This allows the following logic to fire,
  //   and prevent an additional ApplyBoundaryAction for such turnaround steps.
  G4LatticePhysical* volLattice =
    G4LatticeManager::GetLatticeManager()->GetLattice(prePV);
  G4cout << "In G4CMPBoundaryUtils: PreStepPoint: " << prePV->GetName() << ", PostStepPoint: " << postPV->GetName() << G4endl;
  G4cout << "In G4CMPBoundaryUtils: VolLattice: " << volLattice << ", procUtils->GetLattice(): " << procUtils->GetLattice() << G4endl;
  if(volLattice != procUtils->GetLattice()){
    G4cout << "\n\n\n REL THE DEFAULT G4CMPBOUNDARYUTILS::GetBoundingVolumes() flag is firing. WHY NOW?" <<G4endl;
    if(buVerboseLevel>1){
      G4cout << procName << ": Track inbound after reflection" << G4endl;
    }
    return false;
  }
  */

  if (buVerboseLevel>1) {
    G4cout <<   "  RELPreStep volume: " << prePV->GetName() << " @ "
	   << aStep.GetPreStepPoint()->GetPosition()
	   << "\n PostStep volume: " << (postPV?postPV->GetName():"OutOfWorld")
	   << " @ " << aStep.GetPostStepPoint()->GetPosition()
       << "\n PreStep momentum direction"<< aStep.GetPreStepPoint()->GetMomentumDirection()
       << "\n PostStep momentum direction"<< aStep.GetPostStepPoint()->GetMomentumDirection()
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
  if (G4CMP::IsBogoliubovQP(pd)) {
    matTable = surfProp->GetBogoliubovQPMaterialPropertiesTablePointer();
    electrode = surfProp->GetBogoliubovQPElectrode();
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
  GetBoundingVolumes(aStep); //REL uhhhh does this... do anything?
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
//  G4cout << "postIn (in CheckStepBoundary): " << postIn << G4endl;

  
  if (buVerboseLevel>2) {
    G4cout << "\n Is postStep location on surface of preStep Volume? "
	   << (postIn==kOutside ? "outside" :
	       postIn==kInside  ? "inside" :
	       postIn==kSurface ? "surface" : "INVALID") << G4endl;
  }

  // If post-step position not proper surface point, compute intersection
  if (postIn != kSurface) {
    G4ThreeVector along = (postPos-prePos).unit();	// Trajectory direction
    surfacePoint = prePos + preSolid->DistanceToOut(prePos,along)*along;

    // Double check calculation -- point "must" now be on surface!
//    G4cout << "REL preSolid->Inside(surfacePoint): " << preSolid->Inside(surfacePoint) << G4endl;
//    G4cout << "REL kInside: " << kInside << ", kSurface: " << kSurface << ", kOutside: " << kOutside << G4endl;
    if (preSolid->Inside(surfacePoint) != kSurface) {
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

  bool trackedSCResponse = true;
  if (!matTable) {
    DoSimpleKill(aTrack, aStep, aParticleChange);
  } else if (electrode && electrode->IsNearElectrode(aStep) && !trackedSCResponse) {
    electrode->AbsorbAtElectrode(aTrack, aStep, aParticleChange);
  } else if (AbsorbTrack(aTrack, aStep)) {
    G4cout << "DoingAbsorption in ApplyBoundaryAction" << G4endl;
    DoAbsorption(aTrack, aStep, aParticleChange);
  } else if (MaximumReflections(aTrack)) {
    G4cout << "Maximum reflections reached in ApplyBoundaryAction" << G4endl;
    DoSimpleKill(aTrack, aStep, aParticleChange);
  } else if (ReflectTrack(aTrack, aStep)) {
    G4cout << "Doing reflection in ApplyBoundaryAction" << G4endl;
    IncrementReflectionCount(aTrack);
    DoReflection(aTrack, aStep, aParticleChange);
  } else {
    G4cout << "Doing transmission in ApplyBoundaryAction" << G4endl;
    DoTransmission(aTrack, aStep, aParticleChange);
  }
}


//Dedicated function for doing this. I don't think this should exist in the "check" functions, since
//it will run even if there is no reflection at a surface (i.e. if there is transmission). Putting
//in own function so that we don't have to spread it among a bunch of derived class functions.
void G4CMPBoundaryUtils::IncrementReflectionCount(const G4Track& aTrack)
{
    auto trackInfo = G4CMP::GetTrackInfo<G4CMPVTrackInfo>(aTrack);
    trackInfo->IncrementReflectionCount();
}

// Default conditions for absorption or reflection

G4bool G4CMPBoundaryUtils::AbsorbTrack(const G4Track&, const G4Step&) const {
  G4double absProb = GetMaterialProperty("absProb");
    
  if (buVerboseLevel>2) G4cout << " AbsorbTrack: absProb " << absProb << G4endl;

  return (G4UniformRand() <= absProb);
}

G4bool G4CMPBoundaryUtils::ReflectTrack(const G4Track& aTrack, const G4Step&) const {
  G4double reflProb = GetMaterialProperty("reflProb");
    
  if (buVerboseLevel>2) G4cout << " ReflectTrack: reflProb " << reflProb << G4endl;

  return (G4UniformRand() <= reflProb);
}

G4bool G4CMPBoundaryUtils::MaximumReflections(const G4Track& aTrack) const {
  auto trackInfo = G4CMP::GetTrackInfo<G4CMPVTrackInfo>(aTrack);
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

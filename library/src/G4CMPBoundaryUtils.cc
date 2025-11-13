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
     : G4CMP::IsPhonon(pd) ? G4CMPConfigManager::GetMaxPhononBounces()
     : G4CMP::IsQP(pd) ? G4CMPConfigManager::GetMaxQPBounces() : -1);

  if (buVerboseLevel>1) {
    G4cout << procName << "::IsGoodBoundary maxRefl " << maximumReflections
	   << G4endl;
  }

  return (IsBounaryStep(aStep) &&
	  GetBoundingVolumes(aStep) &&
	  GetSurfaceProperty(aStep));
}

G4bool G4CMPBoundaryUtils::IsBounaryStep(const G4Step& aStep) {

  //Debugging
  if (buVerboseLevel > 5) {
    G4cout << "---------- G4CMPBoundaryUtils::IsBounaryStep ----------" << G4endl;
    G4cout << "IBS Function Point A | The step status is " << aStep.GetPostStepPoint()->GetStepStatus() << G4endl;
  }
  
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

  //Debugging
  if (buVerboseLevel > 5) {
    G4cout << "---------- G4CMPBoundaryUtils::GetBoundingVolumes ----------" << G4endl;
  }

  
  prePV = aStep.GetPreStepPoint()->GetPhysicalVolume();
  postPV = aStep.GetPostStepPoint()->GetPhysicalVolume();

  if (prePV == postPV) {
    if (buVerboseLevel) {
      G4cerr << procName << " ERROR: fGeomBoundary status set, but"
	     << " pre- and post-step volumes are identical!" << G4endl;
    }

    //Debugging
    if (buVerboseLevel > 5) {
      G4cout << "GBV Function Point A | prePV and postPV are the same." << G4endl;
    }
    return false;
  }

  //Note that since now the GetMFP function has called a lattice update, procUtils->GetLattice() should just be the "current" lattice we're in.
  //Debugging
  if (buVerboseLevel > 5) {
    G4cout << "GBV Function Point B | procUtils->GetLattice(): " << procUtils->GetLattice() << ", prePVlattice: " << G4LatticeManager::GetLatticeManager()->GetLattice(prePV) << ", volLattice (postPVLattice): " << G4LatticeManager::GetLatticeManager()->GetLattice(postPV) << G4endl;
  }
  G4LatticePhysical* volLattice = G4LatticeManager::GetLatticeManager()->GetLattice(postPV);

  //First scenario: we DON'T have a lattice in the pre-step point. This happens at least when we're on the boundary of a lattice and    
  //the world. Here the phonon "turns around" in the world volume, in a step that is of negligible length. This negligible length is   
  //typically around 1E-15, so to account for these safely we use a tolerance of 1E-13, which is well below the physics scales       
  //relevant in these kinds of sims.
  double stepLengthTolerance = 1E-13 * CLHEP::m;
  if (G4LatticeManager::GetLatticeManager()->GetLattice(prePV) == 0 ) {
    if (buVerboseLevel > 5) {
      G4cout << "GBV Function Point C | prePV lattice is zero." << G4endl;
    }
    
    //First: if the current (i.e. pre-step, procUtils->GetLattice()) lattice is the same as post-step volume lattice.
    //This occurs, if, for example, the current volume is World/vacuum, and a lattice changeover/update failed in the MFP step because
    //there is no lattice to update to. Hence, the "current" procUtils->GetLattice() lattice is really the prePV one from the *previous* step.
    if (volLattice == procUtils->GetLattice()) {
      if (buVerboseLevel > 5 ) {
	G4cout << "GBV Function Point D | Current (procUtils) Lattice is equal to post-PV lattice." << G4endl;
      }
      
      //If the step length is below tolerance, we need to return false so we don't try to "double-count" the boundary action.
      //The small step sizes occur when the phonon "turns around" on a boundary with a volume that doesn't have a lattice.
      if (aStep.GetStepLength() <= stepLengthTolerance) {
	if (buVerboseLevel > 5) {
	  G4cout << "GBV Function Point E | Step length, " << aStep.GetStepLength()*1.0e9 << " (mult x 1e9) is below step length tolerance." << G4endl;
	}
        return false;
      }
      //If the step length is long. I think this can only happen if we're somehow having a long track in a non-lattice volume. For now
      //we'll keep the control block here just in case (and make it return what would be returned in its absence anyway).
      else{
	if (buVerboseLevel > 5) {
	  G4cout << "GBV Function Point F | Step length is above step length tolerance." << G4endl;
	}
        return true; //TBD, but set to be consistent with older version
      }
    }
    //If the current lattice is not equal to the post-step lattice. From a first glance at phonon dynamics I don't think this happens,
    //but we'll keep the control block here so that if we see it does happen, we can make the call then.
    else{
      if (buVerboseLevel > 5) {
	G4cout << "GBV Function Point G | Current Lattice is not equal to post-PV lattice." << G4endl;
      }
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
  if (G4LatticeManager::GetLatticeManager()->GetLattice(prePV) != 0 ) {
    if (buVerboseLevel > 5) {
      G4cout << "GBV Function Point H | Current lattice is not zero." << G4endl;
    }
    
    //If the current (i.e. pre-step, procUtils->GetLattice()) lattice is different from the post-step volume lattice, then this is where
    //some important logic must happen.
    if (volLattice != procUtils->GetLattice()) {
      if (buVerboseLevel > 5) {
	G4cout << "GBV Function Point I | Current Lattice ("<< procUtils->GetLattice() << ") is not equal to post-PV lattice (" << volLattice << ")." << G4endl;
      }

      //If the step length is tiny, this means that we're in a "turnaround" step characteristic of reflection back into the lattice
      //that we approached from. Here we return false so that we don't try to run another boundary process action for this turnaround
      //step -- that has already been done
      if (aStep.GetStepLength() <= stepLengthTolerance) {
	if (buVerboseLevel > 5) {	  
	  G4cout << "GBV Function Point J | Step length, " << aStep.GetStepLength()*1.0e9 << " (mult x 1e9), is below step length tolerance." << G4endl;
	}
        return false;
      }
      //Otherwise, the track is actually moving through a volume before it hits this surface, and it needs to run the logic to see if
      //it reflects, transmits, etc. Need to return true so that logic can run.
      else{
	if (buVerboseLevel > 5) {
	  G4cout << "GBV Function Point K | Step length, " << aStep.GetStepLength()*1.0e9 << " (mult x 1e9), is above step length tolerance." << G4endl;
	}
        return true;
      }
    }
    //If the current lattice is NOT different from the post-step lattice, this should have already been covered at the beginning of this
    //function -- throw an error here.
    else{
      if (buVerboseLevel > 5) {
	G4cout << "GBV Function Point L | Current lattice is not different from postPV lattice. This should have been covered already -- figure out why it hasn't.\n\n\n" << G4endl;
      }
      return false; //TBD
    }
  }

  if (buVerboseLevel>1) {
    G4cout <<   "  PreStep volume: " << prePV->GetName() << " @ "
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
  if (G4CMP::IsQP(pd)) {
    matTable = surfProp->GetQPMaterialPropertiesTablePointer();
    electrode = surfProp->GetQPElectrode();
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
// "surfacePoint" returns post-step position, or computed surface point. Extra logic being added
// in trackedFilmResponse is the folowing (and is hopefully sufficiently general...):
// 1. Check pre-PV and post-PV volumes to see if it's on either volume's surface.
// 2. If it is, then no action is needed -- can return true.
// 3. If it's not, then figure out which surface is closer to the current point (using DistanceToIn/DistanceToOut),
//    and then do the calculation that is currently done there.

G4bool G4CMPBoundaryUtils::CheckStepBoundary(const G4Step& aStep,
					     G4ThreeVector& surfacePoint) {
  //Debugging
  if (buVerboseLevel > 5) {
    G4cout << "---------- G4CMPBoundaryUtils::CheckStepBoundary ----------" << G4endl;
  }

  
  G4StepPoint* preP = aStep.GetPreStepPoint();
  G4StepPoint* postP = aStep.GetPostStepPoint();
  GetBoundingVolumes(aStep); 
  surfacePoint = postP->GetPosition();		// Correct if valid boundary

  
  // Get pre- and post-step positions in pre-step volume coordinates
  G4VSolid* preSolid = prePV->GetLogicalVolume()->GetSolid();
  G4ThreeVector prePos = preP->GetPosition();
  G4ThreeVector postPos = surfacePoint;

  //Debugging
  if (buVerboseLevel > 5) {
    G4cout << "CSB Function Point A | postPos, pre-rotation = " << postPos << G4endl;
    G4cout << "CSB Function Point A | surfacePoint, pre-rotation = " << surfacePoint << G4endl;
  }
  G4CMP::RotateToLocalPosition(preP->GetTouchable(), prePos);
  G4CMP::RotateToLocalPosition(preP->GetTouchable(), postPos);
  if (buVerboseLevel > 5) {
    G4cout << "CSB Function Point B | postPos, post-rotation = " << postPos << G4endl;
    G4cout << "CSB Function Point B | surfacePoint, post-rotation = " << surfacePoint << G4endl;
  }
  
  //Get pre- and post-step positions in post-step volume coordinates (for good measure)
  G4VSolid* postSolid = postPV->GetLogicalVolume()->GetSolid();
  G4ThreeVector prePos_postPV = preP->GetPosition();
  G4ThreeVector postPos_postPV = surfacePoint;

  //Debugging
  if (buVerboseLevel > 5) {
    G4cout << "CSB Function Point C | postPos_postPV, pre-rotation = " << postPos_postPV << G4endl;
    G4cout << "CSB Function Point C | surfacePoint, pre-rotation = " << surfacePoint << G4endl;
  }
  G4CMP::RotateToLocalPosition(postP->GetTouchable(), prePos_postPV);
  G4CMP::RotateToLocalPosition(postP->GetTouchable(), postPos_postPV);
  if (buVerboseLevel > 5) {
    G4cout << "CSB Function Point D | postPos_postPV, post-rotation = " << postPos_postPV << G4endl;
    G4cout << "CSB Function Point D | surfacePoint, post-rotation = " << surfacePoint << G4endl;
  }
  
  if (buVerboseLevel>2) {
    G4cout << "CheckStepBoundary: in prePV (" << prePV->GetName() << ") frame"
	   << "\n  preStep @ " << prePos << "\n postStep @ " << postPos
	   << G4endl;
    G4cout << "CheckStepBoundary: in postPV (" << postPV->GetName() << ") frame"
	   << "\n  preStep @ " << prePos_postPV << "\n postStep @ " << postPos_postPV
	   << G4endl;
  }

  // Verify that post-step position is on surface of volume. Check both the mother and daughter
  // volumes, since we may have nested geometry.
  EInside postIn = preSolid->Inside(postPos);
  EInside postIn_postPV = postSolid->Inside(postPos_postPV);
  
  if (buVerboseLevel>2) {
    G4cout << "\n Is postStep location on surface of preStep Volume? "
	   << (postIn==kOutside ? "outside" :
	       postIn==kInside  ? "inside" :
	       postIn==kSurface ? "surface" : "INVALID") << G4endl;
    G4cout << "\n Is postStep location on surface of postStepVolume? "
	   << (postIn_postPV==kOutside ? "outside" :
	       postIn_postPV==kInside  ? "inside" :
	       postIn_postPV==kSurface ? "surface" : "INVALID") << G4endl;
  }

  // If post-step position is on the post-PV boundary OR if post-step position is on the pre-PV boundary,
  // then we can return true.
  if (postIn == kSurface || postIn_postPV == kSurface) {
    return true;
  }
  
  //Otherwise, we need a bit of logic to handle the adjustment, since it's not obvious which boundary it should be targeting.
  //REL: This is really nasty. Can we do the following logicking in a way that's more transparent? Basically, need to:
  //1. Identify which volume (pre/post-PV) is internal to the other
  //2. Identify where the post-step point is relative to the volumes
  //3. Identify which volume (pre/post-PV) is closer to the post-step point
  //Unfortunately this gives a nested loop with 2^3 potential options, which is gross and error-prone. Is there a nicer way to
  //get this information so we can do the adjustments if we're not on a surface?
  
  
  //---------------------------------------------------------------------------------------------------------------------
  //If we're outside the pre-step volume, check whether we're inside or outside the post-step volume
  if (postIn == kOutside) {

    //-------------------------------
    //If we're also outside the post-step volume, uhhhhh... where are we?
    if (postIn_postPV == kOutside) {
      G4Exception((procName+"::CheckStepBoundary").c_str(),"Boundary00X",JustWarning,
		  ("Post-step point is somehow outside both the pre-step volume, " +
		   prePV->GetName() +
		   " and the post-step volume, " +
		   postPV->GetName()).c_str());

	/////
      /*  REL this stuff is from NT I think. Chat with NT/MK about whether this is 100% necessary given TFR's upgrade
	  But otherwise I think I will keep mine because it makes more general assumptions about phonons impinging on
	  the interface between two crystals
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
      */
      return false;
    }

    //-------------------------------
    //Otherwise, we're outside pre-step volume and inside the post-step volume -- this is sensible if:
    //1. The pre-step volume is the daughter volume, the post-step volume is the mother volume, and the point is between the two boundaries...
    //2. ...or if the volumes are not nested (i.e. neither is a daughter of the other)
    else if (postIn_postPV == kInside) {

      //Here, compare how far inside the post-step volume we are with how far outside the pre-step volume we are.
      //Not using "along" vector here for simplicity
      G4double pre_distToIn = preSolid->DistanceToIn(postPos);
      G4double post_distToOut = postSolid->DistanceToOut(postPos_postPV);

      //------------
      //We're closer to the pre-step volume's surface than the post-step volume's surface.
      if (fabs(post_distToOut) > fabs(pre_distToIn)) {

	//Put point onto the pre-step volume surface, now using the distanceToIn with the right direction.
	//The post-step point is outside the pre-PV and "along" points out. Need to subtract off outward vector (hence minus sign)
	G4ThreeVector along = (postPos-prePos).unit(); // Trajectory direction	
	surfacePoint = postPos - fabs(preSolid->DistanceToIn(postPos,along)) * along; 

	//Check that the surface point is good (now that we've modified it, it's in the local coordinate system)
	if (preSolid->Inside(surfacePoint) != kSurface) {
	  G4Exception((procName+"::CheckBoundaryPoint").c_str(),
		      "Boundary00X", EventMustBeAborted,
		      "Boundary-limited step cannot find boundary surface point"
		      );
	  return false;
	}
	//Rotate back to global using rotation of the pre-step volume
	G4CMP::RotateToGlobalPosition(preP->GetTouchable(), surfacePoint);
      }

      //------------
      //We're closer to the post-step volume's surface than the pre-step volume's surface, so do the calculation in the post-step volume's coords
      else{

	//Put point onto the post-step volume surface.
	//We're inside the post-PV volume and the along vector points out. Add a positive along vector to get to extrenal boundary
	G4ThreeVector along = (postPos_postPV-prePos_postPV).unit(); // Trajectory direction	
	surfacePoint = postPos_postPV + fabs(postSolid->DistanceToOut(postPos_postPV,along)) * along; 
	
	//Check that the surface point is good (now that we've modified it, it's in the local coordinate system)
	if (postSolid->Inside(surfacePoint) != kSurface) {
	  G4Exception((procName+"::CheckBoundaryPoint").c_str(),
		      "Boundary00X", EventMustBeAborted,
		      "Boundary-limited step cannot find boundary surface point"
		      );
	  return false;
	}

	//Rotate back to global using rotation of the post-step volume
	G4CMP::RotateToGlobalPosition(postP->GetTouchable(), surfacePoint);
      }
	
    }

    //-------------------------------
    //If neither of these get flagged, then we should have something invalid, since boundary cases should have been caught higher up.
    else{
      G4Exception((procName+"::CheckBoundaryPoint").c_str(),
		  "Boundary00X", EventMustBeAborted,
		  "Somehow the post-step point for this step is neither inside, outside, or on the surface of the post-step volume."
		  );
      return false;
    }
  }

  //---------------------------------------------------------------------------------------------------------------------
  //If we're inside the pre-step point's volume, then do a check
  else if (postIn == kInside) {

    //-------------------------------
    //This condition makes sense in the following scenarios:
    //1. The pre-step volume is the mother volume, the post-step volume is the daughter, and the post step point is
    //   inside the former but outside the latter
    //2. The two volumes are not mothers/daughters of each other (i.e. same level of the heirarchy)
    if (postIn_postPV == kOutside) {

      //Here, compare how far inside the post-step volume we are with how far outside the pre-step volume we are.
      G4double pre_distToOut = preSolid->DistanceToOut(postPos);
      G4double post_distToIn = postSolid->DistanceToIn(postPos_postPV);

      //------------
      //If we're closer to the pre-step volume, then put the point on that surface      
      if (fabs(pre_distToOut) < fabs(post_distToIn)) {

	//Put point onto the pre-step volume surface
	//We're inside the pre-step volume and want to get closer to its surface. Along step points inward, but we want to go outward. Need a minus.
	G4ThreeVector along = (postPos-prePos).unit(); // Trajectory direction points in the direction of the step
	surfacePoint = postPos - fabs(preSolid->DistanceToOut(postPos,along)) * along; 

	//Check that the surface point is good (now that we've modified it, it's in the local coordinate system)
	if (preSolid->Inside(surfacePoint) != kSurface) {
	  G4Exception((procName+"::CheckBoundaryPoint").c_str(),
		      "Boundary00X", EventMustBeAborted,
		      "Boundary-limited step cannot find boundary surface point"
		      );
	  return false;
	}

	//Rotate back to global using rotation of the post-step volume
	G4CMP::RotateToGlobalPosition(preP->GetTouchable(), surfacePoint);
	
      }
      //------------
      //If we're closer to the post-step volume, then put the point on that surface.
      else{

	//Put point onto the post-step volume surface
	//We're outside the post-step volume and want to get closer to its surface. Along step points inward. Need a plus.
	G4ThreeVector along = (postPos_postPV-prePos_postPV).unit(); // Trajectory direction	
	surfacePoint = postPos_postPV + fabs(postSolid->DistanceToIn(postPos_postPV,along)) * along; 

	//Check that the surface point is good (now that we've modified it, it's in the local coordinate system)
	if (postSolid->Inside(surfacePoint) != kSurface) {
	  G4Exception((procName+"::CheckBoundaryPoint").c_str(),
		      "Boundary00X", EventMustBeAborted,
		      "Boundary-limited step cannot find boundary surface point"
		      );
	  return false;
	}

	//Rotate back to global using rotation of the post-step volume
	G4CMP::RotateToGlobalPosition(postP->GetTouchable(), surfacePoint);
      }	
    }

    //-------------------------------
    //This condition makes sense (I think?) if the post-step point is within the daughter volume (it's therefore inside both mother and daughter volumes)
    else if (postIn_postPV == kInside) {

      //Since we still don't know which is the daughter, we again have to run a conditional. See which one has a closer boundary (in any direction)
      G4double pre_distToOut = preSolid->DistanceToOut(postPos);
      G4double post_distToOut = postSolid->DistanceToOut(postPos_postPV);

      //Here, the post-PV is internal to the pre-PV
      if (fabs(pre_distToOut) > fabs(post_distToOut)) {

	//Put the point on the post-PV volume, since it's smaller/internal
	//We're inside the post-step volume, which is inside the pre-step volume, and want to land on the post-step surface. Along step points inward. Need
	//a minus to get to post-step surface.
	G4ThreeVector along = (postPos_postPV-prePos_postPV).unit(); // Trajectory direction	
	surfacePoint = postPos_postPV - fabs(postSolid->DistanceToOut(postPos_postPV,along)) * along;

	//Check that the surface point is good (now that we've modified it, it's in the local coordinate system)
	if (postSolid->Inside(surfacePoint) != kSurface) {
	  G4Exception((procName+"::CheckBoundaryPoint").c_str(),
		      "Boundary00X", EventMustBeAborted,
		      "Boundary-limited step cannot find boundary surface point"
		      );
	  return false;
	}

	//Rotate back to global using rotation of the post-step volume
	G4CMP::RotateToGlobalPosition(postP->GetTouchable(), surfacePoint);
      }
      //Here, the pre-PV is internal to the post-PV
      else{
	
	//Put the point on the pre-PV volume, since it's smaller/internal
	//We're inside the pre-step volume, which is inside the post-step volume, and want to land on the pre-step surface. Along step points outward. Need
	//a plus to get to pre-step surface.
	G4ThreeVector along = (postPos-prePos).unit(); // Trajectory direction	
	surfacePoint = postPos + fabs(preSolid->DistanceToOut(postPos,along)) * along; 

	//Check that the surface point is good (now that we've modified it, it's in the local coordinate system)
	if (preSolid->Inside(surfacePoint) != kSurface) {
	  G4Exception((procName+"::CheckBoundaryPoint").c_str(),
		      "Boundary00X", EventMustBeAborted,
		      "Boundary-limited step cannot find boundary surface point"
		      );
	  return false;
	}

	//Rotate back to global using rotation of the post-step volume
	G4CMP::RotateToGlobalPosition(preP->GetTouchable(), surfacePoint);
      }
    }
    //-------------------------------
    //If neither of these get flagged, then we should have something invalid, since boundary cases should have been caught higher up.
    else{
      G4Exception((procName+"::CheckBoundaryPoint").c_str(),
		  "Boundary00X", EventMustBeAborted,
		  "Somehow the post-step point for this step is neither inside, outside, or on the surface of the post-step volume."
		  );
      return false;
    }
  }

  //---------------------------------------------------------------------------------------------------------------------
  //Otherwise, we have an invalid thing?
  else{    
    G4Exception((procName+"::CheckBoundaryPoint").c_str(),
		"Boundary00X", EventMustBeAborted,
		"Somehow the post-step point for this step is neither inside, outside, or on the boundary of the pre-step volume."
		);
    return false;
  }
  
	    
  //If we reach this point, we can return true
  return true;
}



// Implement PostStepDoIt() in a common way; processes should call through

void 
G4CMPBoundaryUtils::ApplyBoundaryAction(const G4Track& aTrack,
					const G4Step& aStep,
					G4ParticleChange& aParticleChange) {

  //Debugging
  if (buVerboseLevel > 5) {
    G4cout << "---------- G4CMPBoundaryUtils::ApplyBoundaryAction ----------" << G4endl;
  }
  
  aParticleChange.Initialize(aTrack);

  if (buVerboseLevel > 5) {
    G4cout << "Track momentum direction: " << aTrack.GetMomentumDirection() << G4endl;
  }

  if (!matTable) {
    if (buVerboseLevel>2) G4cout << "BU::Apply: !matTable" << G4endl;
    DoSimpleKill(aTrack, aStep, aParticleChange);
  } else if (electrode && electrode->IsNearElectrode(aStep) ) {
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
    IncrementReflectionCount(aTrack);
    DoReflection(aTrack, aStep, aParticleChange);
  } else {
    if (buVerboseLevel>2) G4cout << "BU::Apply: Transmission" << G4endl;
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

G4bool G4CMPBoundaryUtils::ReflectTrack(const G4Track& /*aTrack*/, const G4Step&) const {
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
  G4ThreeVector norm = G4CMP::GetSurfaceNormal(aStep,pdir); // Outward normal
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



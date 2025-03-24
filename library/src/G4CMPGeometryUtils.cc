/*********************************************************************** \
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
// File: G4CMPGeometryUtils.cc
//
// Description: Free standing helper functions for geometry based calculations.
//
// 20161107  Rob Agnese
// 20170605  Pass touchable from track, not just local PV
// 20170815  Move AdjustSecondaryPosition to here as ApplySurfaceClearance
// 20170913  Add utility to get electric field at (global) position
// 20170925  Add utility to create touchable at (global) position
// 20190226  Use local instance of G4Navigator to avoid corrupting tracking

#include "G4CMPGeometryUtils.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPGlobalLocalTransformStore.hh"
#include "G4Field.hh"
#include "G4FieldManager.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalVolume.hh"
#include "G4Navigator.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"
#include "G4VTouchable.hh"
#include "G4SystemOfUnits.hh"

G4ThreeVector G4CMP::GetLocalDirection(const G4VTouchable* touch,
                                       const G4ThreeVector& dir) {
  return G4CMPGlobalLocalTransformStore::ToLocal(touch).TransformAxis(dir);
}

G4ThreeVector G4CMP::GetLocalPosition(const G4VTouchable* touch,
                                      const G4ThreeVector& pos) {
  return G4CMPGlobalLocalTransformStore::ToLocal(touch).TransformPoint(pos);
}

G4ThreeVector G4CMP::GetGlobalDirection(const G4VTouchable* touch,
                                        const G4ThreeVector& dir) {
  return G4CMPGlobalLocalTransformStore::ToGlobal(touch).TransformAxis(dir);
}

G4ThreeVector G4CMP::GetGlobalPosition(const G4VTouchable* touch,
                                       const G4ThreeVector& pos) {
  return G4CMPGlobalLocalTransformStore::ToGlobal(touch).TransformPoint(pos);
}

void G4CMP::RotateToLocalDirection(const G4VTouchable* touch,
                                   G4ThreeVector& dir) {
  G4CMPGlobalLocalTransformStore::ToLocal(touch).ApplyAxisTransform(dir);
}

void G4CMP::RotateToLocalPosition(const G4VTouchable* touch,
                                  G4ThreeVector& pos) {
  G4CMPGlobalLocalTransformStore::ToLocal(touch).ApplyPointTransform(pos);
}

void G4CMP::RotateToGlobalDirection(const G4VTouchable* touch,
                                    G4ThreeVector& dir) {
  G4CMPGlobalLocalTransformStore::ToGlobal(touch).ApplyAxisTransform(dir);
}

void G4CMP::RotateToGlobalPosition(const G4VTouchable* touch,
                                   G4ThreeVector& pos) {
  G4CMPGlobalLocalTransformStore::ToGlobal(touch).ApplyPointTransform(pos);
}

/*
//Overloaded version of this for when we only have a single step point (i.e. during MFP calculation at the beginning of a step).
//Since with a single step point we don't necessarily know about multiple volumes, the second argument is a guessed direction that
//can be used to find a nearby volume.
G4ThreeVector G4CMP::GetSurfaceNormal(const G4StepPoint* stepPoint, const G4ThreeVector& guessedDirection )
{
  const G4VTouchable* preTouch = stepPoint->GetTouchable();
  G4VTouchable * guessTouch = CreateTouchableAtPoint(stepPoint->GetPosition()+0.1*nm*guessedDirection); //HARDCODED -- BEWARE  
   
  //For now, throw an exception if the pre-step Volume and the guess volume are the same
  if( preTouch->GetVolume() == guessTouch->GetVolume() ){
    G4ExceptionDescription msg;
    msg << "G4CMP::GetSurfaceNormal()'s preStepVolume seems to be the same as the guess volume, which shouldn't be true. Volumes are: "
	<< preTouch->GetVolume()->GetName() << " and " << guessTouch->GetVolume()->GetName() << G4endl;
    G4Exception("G4CMP::GetSurfaceNormal(stepPoint,guessedDirection)", "Geometry001",
		EventMustBeAborted, msg);
  }

  //Now we use these two volumes with the same logic as before to find the surface normal
  G4VSolid* preSolid = preTouch->GetVolume()->GetLogicalVolume()->GetSolid();
  G4VSolid* guessSolid = guessTouch->GetVolume()->GetLogicalVolume()->GetSolid();

  //Check to see if the position is within spitting distance of either object's boundary. First, understand whether we're
  //inside either of the
  G4ThreeVector pos_prePV = stepPoint->GetPosition();
  G4ThreeVector pos_guessPV = stepPoint->GetPosition();
  G4cout << "pos_prePV1: " << pos_prePV << ", in volume " << preTouch->GetVolume()->GetName() << G4endl;
  G4cout << "pos_guessPV1: " << pos_guessPV << ", in volume " << guessTouch->GetVolume()->GetName() << G4endl;
  
  RotateToLocalPosition(preTouch, pos_prePV);
  RotateToLocalPosition(guessTouch,pos_guessPV);
  EInside postStepInPrePV = preSolid->Inside(pos_prePV);
  EInside postStepInGuessPV = guessSolid->Inside(pos_guessPV);
  G4cout << "pos_prePV2: " << pos_prePV << G4endl;
  G4cout << "pos_guessPV2: " << pos_guessPV << G4endl;

  
  G4cout << "PostStepInPrePV: " << postStepInPrePV << ", postStepInGuessPV: " << postStepInGuessPV << ", fabs(preSolid->DistanceToOut(pos_prePV)):" << fabs(preSolid->DistanceToOut(pos_prePV)) << ", fabs(preSolid->DistanceToIn(pos_prePV)): " << fabs(preSolid->DistanceToIn(pos_prePV)) << ", fabs(guessSolid->DistanceToOut(pos_guessPV)): " << fabs(guessSolid->DistanceToOut(pos_guessPV)) << ", fabs(guessSolid->DistanceToIn(pos_guessPV)): " << fabs(guessSolid->DistanceToIn(pos_guessPV)) << G4endl;

  
  //Now, we have some logic. We don't want as much logic as in G4CMPBoundaryUtils::CheckBoundarySurface(), especially
  //when it comes to what happens when we don't have a point on a bonafide surface -- that code should be run
  //first to confirm that we're indeed on a boundary surface. But we do need something to tell us which surface's normal to reflect over.
  double tolerance = 1.0e-11 * mm; //Need to not hardcode this -- find a better way to implement. 
  if( postStepInPrePV == kSurface ){
    //Reflect over the pre-PV normal
    G4ThreeVector preSolidNorm = preSolid->SurfaceNormal(pos_prePV);
    G4cout << "REL in G4CMPGeometryUtils::GetSurfaceNormal() --- returning pre-PV surface norm: " << preSolidNorm << G4endl;
    RotateToGlobalDirection(preTouch,preSolidNorm);
    G4cout << "REL in G4CMPGeometryUtils::GetSurfaceNormal() --- after rotation, this surface norm is: " << preSolidNorm << G4endl;
    delete guessTouch;
    return preSolidNorm;
  }
  else if( postStepInGuessPV == kSurface ){
    //Reflect over the guess-PV normal
    G4ThreeVector guessSolidNorm = guessSolid->SurfaceNormal(pos_guessPV);
    G4cout << "REL in G4CMPGeometryUtils::GetSurfaceNormal() --- returning guess-PV surface norm: " << guessSolidNorm << G4endl;
    RotateToGlobalDirection(guessTouch,guessSolidNorm);
    G4cout << "REL in G4CMPGeometryUtils::GetSurfaceNormal() --- after rotation, this surface norm is: " << guessSolidNorm << G4endl;
    return guessSolidNorm;
  }
  //If we are very near a surface and within tolerance, still okay -- this is a pre-PV reflection
  else if( (fabs(preSolid->DistanceToOut(pos_prePV)) > 0 && fabs(preSolid->DistanceToOut(pos_prePV)) < tolerance ) ||
	   (fabs(preSolid->DistanceToIn(pos_prePV)) > 0 && fabs(preSolid->DistanceToIn(pos_prePV)) < tolerance ) ){
    G4ThreeVector preSolidNorm = preSolid->SurfaceNormal(pos_prePV);
    G4cout << "REL in G4CMPGeometryUtils::GetSurfaceNormal() --- returning pre-PV surface norm: " << preSolidNorm << G4endl;
    RotateToGlobalDirection(preTouch,preSolidNorm);
    G4cout << "REL in G4CMPGeometryUtils::GetSurfaceNormal() --- after rotation, this surface norm is: " << preSolidNorm << G4endl;
    delete guessTouch;
    return preSolidNorm;
  }
  //If we are very near a surface and within tolerance, still okay -- this is for guess-PV reflection
  else if( (fabs(guessSolid->DistanceToOut(pos_guessPV)) > 0 && fabs(guessSolid->DistanceToOut(pos_guessPV)) < tolerance ) ||
	   (fabs(guessSolid->DistanceToIn(pos_guessPV)) > 0 && fabs(guessSolid->DistanceToIn(pos_guessPV)) < tolerance ) ){
    G4ThreeVector guessSolidNorm = guessSolid->SurfaceNormal(pos_guessPV);
    G4cout << "REL in G4CMPGeometryUtils::GetSurfaceNormal() --- returning guess-PV surface norm: " << guessSolidNorm << G4endl;
    RotateToGlobalDirection(guessTouch,guessSolidNorm);
    G4cout << "REL in G4CMPGeometryUtils::GetSurfaceNormal() --- after rotation, this surface norm is: " << guessSolidNorm << G4endl;
    delete guessTouch;
    return guessSolidNorm;
  }
  //Otherwise, we're not on (or within tolerance of) a surface.
  else{
    
    //Throw an error -- shouldn't ever be here since it means we're not on a surface
    G4ExceptionDescription msg;
    msg << "G4CMP::GetSurfaceNormal() seems to think we're not on a surface. Volumes are: "
	<< preSolid->GetName() << " and " << guessSolid->GetName() << G4endl;
    G4Exception("G4CMP::GetSurfaceNormal()", "Geometry00X",
		EventMustBeAborted, msg);
    G4ThreeVector dummy(0,0,0);
    delete guessTouch;
    return dummy;    
  }

  
}
*/




//Nicer, concise version of Get2DSafety. Arguments are:
//1. The touchable of the mother volume that the point is currently in
//2. The position of the point whose safety is to be computed
//3. The momentum direction of the last point (used for identifying what directions are "in-plane")
//4. A bool saying whether we're trying compute a safety from a boundary we're currently on. REL do we actually need this, or can we not calculate?
//5. If we are calculating from a boundary, the surface norm at that point.
//6. If we are computing a constrained safety (to get a QP unstuck from a corner), an outgoing tangent vector at the corner
//7. If we are computing a constrained safety (to get a QP unstuck from a corner), the other outgoing tangent vector at the corner
G4double G4CMP::Get2DSafety(const G4VTouchable* motherTouch,
			    G4ThreeVector pos,
			    G4ThreeVector momDir,
			    bool safetyFromABoundary,
			    G4ThreeVector surfaceNorm,
			    G4ThreeVector tangVect1,
			    G4ThreeVector tangVect2)
{
  G4int verboseLevel = G4CMPConfigManager::GetVerboseLevel();
  if( verboseLevel > 5 ){
    G4cout << "---------- G4CMPGeometryUtils::Get2DSafety ----------" << G4endl;
  }



  
  //Pseudocode
  //0. Define output variables
  G4double overallSafety = DBL_MAX;

  //Get the mother volume information
  G4VPhysicalVolume* motherPhys = motherTouch->GetVolume();
  G4LogicalVolume* motherLog = motherPhys->GetLogicalVolume();
  G4VSolid * motherSolid = motherLog->GetSolid();

  //Rotate the position and the momentum direction to the touchable's coordinates. If a surfacenorm, tangVect1, and tangVect2 are not
  //supplied, then these just rotate from the 0 vector to the 0 vector.
  RotateToLocalPosition(motherTouch, pos);
  RotateToLocalDirection(motherTouch, momDir);
  RotateToLocalDirection(motherTouch,surfaceNorm);
  RotateToLocalDirection(motherTouch,tangVect1);
  RotateToLocalDirection(motherTouch,tangVect2);


  //1. First, get the shortest distance to the mother volume that we're in. ("DistanceToOut")
  G4double motherSafety = Compute2DSafetyInMotherVolume(motherSolid, pos, safetyFromABoundary, surfaceNorm, tangVect1, tangVect2);
  if( motherSafety < overallSafety ){
    overallSafety = motherSafety;

    //Debugging
    if( verboseLevel > 5 ){
      G4cout << "G2DS Function Point A | overall safety during mother safety check: " << overallSafety << G4endl;
    }
  }

  //2. Next, loop through the daughter volumes in this mother volume and compute the distances to those.
  //When we're stuck, we don't actually use the tangVector1 and tangVector2 in the math, because of something that I think will be a reasonable
  //first approximation: if you're in a corner, the chances of a daughter boundary being on the same scale of how close you are to that
  //corner should be small. Some edge cases exist, where you can form a corner from a mother and a daugher boundary, in which case you
  //may land on the daughter, re-recognize that you're stuck, and then get re-ejected safely far from the corner. So at least in reasonably simple
  //geometries, this plus our mother safety should be fine. Need to rejigger this code to integrate it into the Get2D safety, since there are lots
  //of similarities/overlaps. (REL)
  for( int iD = 0; iD < motherLog->GetNoDaughters(); ++iD ){
    G4double daughterSafety = Compute2DSafetyToDaughterVolume(pos,momDir,motherLog,safetyFromABoundary,iD,surfaceNorm,tangVect1,tangVect2);
    if( daughterSafety < overallSafety ) overallSafety = daughterSafety;

    //Debugging
    if( verboseLevel > 5 ){
      G4cout << "G2DS Function Point B | Overall safety during daughter safety check " << iD << ": " << overallSafety << G4endl;
    }
  }
  
  //Return the safety
  return overallSafety; 
  
}

//Looking "inward" within a mother volume to its daughter volumes to identify safeties. Generally, don't want to use this on its
//own. Should only be called from Get2DSafety, not by separate classes.
G4double G4CMP::Compute2DSafetyToDaughterVolume(const G4ThreeVector & pos, const G4ThreeVector & momDir, G4LogicalVolume * motherLog, bool safetyFromABoundary, G4int daughterID, G4ThreeVector surfaceNorm, G4ThreeVector tangVect1, G4ThreeVector tangVect2 ){

  G4int verboseLevel = G4CMPConfigManager::GetVerboseLevel();
  if( verboseLevel > 5 ){
    G4cout << "---------- G4CMPGeometryUtils::Compute2DSafetyToDaughterVolume ----------" << G4endl;
  }

  
  //Establish output variable;
  G4double safety = DBL_MAX;
  
  //Get the physical volume of the daughter from the mother logical and the iD
  G4VPhysicalVolume * volDaughterPhys = motherLog->GetDaughter(daughterID);
  
  //This stuff is ripped verbatim from G4NormalNavigation. Rotates the point to be in the
  //standard reference frame of our daughter's G4solid
  G4AffineTransform sampleTf(volDaughterPhys->GetRotation(),
			     volDaughterPhys->GetTranslation());
  
  sampleTf.Invert();
  const G4ThreeVector samplePoint = sampleTf.TransformPoint(pos);
  const G4ThreeVector sampleAxis = sampleTf.TransformAxis(momDir);
  const G4ThreeVector rotatedSurfaceNorm = sampleTf.TransformAxis(surfaceNorm);
  const G4ThreeVector rotatedTangVect1 = sampleTf.TransformAxis(tangVect1);
  const G4ThreeVector rotatedTangVect2 = sampleTf.TransformAxis(tangVect2);
  const G4VSolid * volDaughterSolid = volDaughterPhys->GetLogicalVolume()->GetSolid();

  //Check to make sure we're infact outside the daughter volume
  if( volDaughterSolid->Inside(samplePoint) == kInside ){
    G4ExceptionDescription msg;
    msg << "G4CMP::Compute2DSafetyToDaughterVolume seems to think we're already inside the daughter volume." << G4endl;
    G4Exception("G4CMP::Compute2DSafetyToDaughterVolume()", "Geometry00X",FatalException, msg);
  }

  //Compute the safety. Since we should be *outside* the daughter volume laterally, it means that the distToIn
  //is in fact the distance in XY, assuming the distanceToIn is computed correctly. This is not exactly true for
  //a few solids, which will require upgrades. (See G4Box)
  const G4double volDaughterSafety = volDaughterSolid->DistanceToIn(samplePoint);
  
  //More complicated logic, which takes into account the effect of being on a boundary of a surface with
  //concavities.
  if( safetyFromABoundary ){

    //Debugging
    if( verboseLevel > 5 ){
      G4cout << "G2DSTDV Function Point A | SafetyFromABoundary is triggering, with volDaughterSafety: " << volDaughterSafety << G4endl;
    }



    //Determine if the current boundary belongs to this daughter
    G4bool currentBoundaryBelongsToThisDaughter = false;
    G4double boundaryTolerance = 1e-12; //REL HARDCODED NEED TO FIX -- maybe use kSurface here instead?
    if( volDaughterSafety < boundaryTolerance ){
      currentBoundaryBelongsToThisDaughter = true;

      //Debugging
      if( verboseLevel > 5 ){
	G4cout << "G2DSTDV Function Point B | Changing currentBoundaryBelongsToThisDaughter to true." << G4endl;
      }
    }
    
    //If it doesn't, then we can just proceed as normal -- I think there aren't odd edge cases here
    if( !currentBoundaryBelongsToThisDaughter ){
      if( volDaughterSafety < safety ) safety = volDaughterSafety;
    }    
    //If it does, we have do do a sweep so we can see if there are places on this daughter volume
    //that are closer to the given point than the safeties to other daughter volumes (and the mother)
    else{
      G4double safetyToThisDaughterBoundary = Compute2DSafetyFromABoundary(volDaughterSolid,samplePoint,rotatedSurfaceNorm,rotatedTangVect1,rotatedTangVect2,false);
      if( safetyToThisDaughterBoundary < safety ) safety = safetyToThisDaughterBoundary;
    }
  }
  //If we're not on a boundary we get easy logic
  else{
    if( volDaughterSafety < safety ) safety = volDaughterSafety;
  }  

  //Debugging
  if( verboseLevel > 5 ){
    G4cout << "G2DSTDV Function Point C | we are looking at a daughter volume of: " << volDaughterPhys->GetName() << " which has a distToIn of " << safety << G4endl;
  }
  return safety;
}



//Looking "outward" from the volume we're in, compute the 2D safety in XY. We do this with a loop, but this will
//need to be replaced with something that's a bit less geometry-agnostic in order to be efficient. This is the place
//to start that further optimization. Generally, don't want to use this on its own. Should only be called from
//Get2DSafety, not by separate classes.
G4double G4CMP::Compute2DSafetyInMotherVolume(G4VSolid * motherSolid,
					      G4ThreeVector pos,
					      bool safetyFromABoundary,
					      G4ThreeVector surfaceNorm,
					      G4ThreeVector tangVect1,
					      G4ThreeVector tangVect2)
{
  G4double motherSafety = DBL_MAX;
  
  //Check to make sure we're in fact inside the mother volume
  if( motherSolid->Inside(pos) == kOutside ){
    G4ExceptionDescription msg;
    msg << "G4CMP::Compute2DSafetyInMotherVolume seems to think we're outside the mother volume." << G4endl;
    G4Exception("G4CMP::Compute2DSafetyInMotherVolume()", "Geometry00X",FatalException, msg);
  }


  //Note that we need to confirm that the "safetyFromABoundary" implies that it's the *mother's* boundary, and
  //not a daughter boundary. As a result, we need to ensure that kSurface is true for the mother here.
  //If we're computing from a boundary, want to standardize the spammed rays with respect to the surface norm
  if( safetyFromABoundary && motherSolid->Inside(pos) == kSurface ){
    G4double motherSafetyFromABoundary = Compute2DSafetyFromABoundary(motherSolid,pos,surfaceNorm,tangVect1,tangVect2,true);
    motherSafety = motherSafetyFromABoundary;
  }
  //Otherwise, doesn't matter -- this is freeform
  else{
    G4double motherSafetyFromTheBulk = Compute2DMotherSafetyFromtheBulk(motherSolid,pos);
    motherSafety = motherSafetyFromTheBulk;
  }

  return motherSafety;
}

//Compute the 2D safety to the mother from a point in the bulk (i.e. not on a boundary)
G4double G4CMP::Compute2DMotherSafetyFromtheBulk(const G4VSolid * motherSolid,G4ThreeVector pos)
{
  G4int verboseLevel = G4CMPConfigManager::GetVerboseLevel();
  if( verboseLevel > 5 ){
    G4cout << "---------- G4CMPGeometryUtils::Compute2DMotherSafetyFromtheBulk ----------" << G4endl;
  }
  
  //Spam a set of DistToOuts in different directions. This part slows things down substantially and should be used sparingly. This part
  //should be replaced with geometry math
  double motherSafety = DBL_MAX;
  G4ThreeVector theDir;

  G4int nV = G4CMPConfigManager::GetSafetyNSweep2D();
  clock_t timestampStart, timestampEnd;
  timestampStart = clock();
  for( G4int iV = 0; iV < nV; ++iV ){
    G4double phi = 2*CLHEP::pi*((double)iV/(double)nV);
    G4double x = cos(phi);
    G4double y = sin(phi);
    theDir.setX(x);
    theDir.setY(y);
    theDir.setZ(0);
    G4double distToOut = motherSolid->DistanceToOut(pos,theDir);
    if( distToOut < motherSafety ) motherSafety = distToOut;
  }
  timestampEnd = clock();

  //Debugging
  if( verboseLevel > 5 ){
    G4cout << "C2DMSFTB Function Point A | Time elapsed during 2D DistToOutLoop: " << double(timestampEnd-timestampStart)/double(CLOCKS_PER_SEC) << " seconds" << G4endl;
    G4cout << "C2DMSFTB Function Point A | clocks_per_sec: " << CLOCKS_PER_SEC << G4endl;
    G4cout << "C2DMSFTB Function Point A | In Compute2DMotherSafetyFromTheBulk, looking at a mother safety of: " << motherSafety << G4endl;
  }
  return motherSafety;
}


//Spam a set of DistToOuts in different directions. This part slows things down substantially and should be used sparingly. This part
//should be replaced with geometry math
G4double G4CMP::Compute2DSafetyFromABoundary(const G4VSolid * theVolSolid, G4ThreeVector pos, G4ThreeVector surfaceNorm, G4ThreeVector tangVect1, G4ThreeVector tangVect2,bool volIsMother)
{
  G4int verboseLevel = G4CMPConfigManager::GetVerboseLevel();
  if( verboseLevel > 5 ){
    G4cout << "---------- G4CMPGeometryUtils::Compute2DSafetyFromABoundary ----------" << G4endl;
  }
  
  double the2DSafety = DBL_MAX;
  G4ThreeVector theDir;
  clock_t timestampStart, timestampEnd;
  timestampStart = clock();
  G4double smallTolerance = 1e-10; //Just in case we're hitting floating point errors.

  //Want to generate our points such that we start from the norm and work outwards
  G4int nV = G4CMPConfigManager::GetSafetyNSweep2D();
  G4double dotProductThreshold_Norm = ComputeDotProductThreshold_Norm(nV);
  G4double nVHalfCircle = nV/2;
  
  //Rotate in positive angular direction
  for( G4int iV = 0; iV < nVHalfCircle; ++iV ){
    G4double deltaPhi = 1.0*CLHEP::pi*((double)iV/(double)nVHalfCircle);
    theDir = surfaceNorm;
    theDir.rotateZ(deltaPhi); //Needs to be plane-agnostic at some point REL

    //Check to make sure that the scan vector is sufficiently away from the surface tangent so we don't end up with
    //ridiculously crazy small steps'
    //G4double epsilonDotProductForNorm = 0.0896393089; //0.07; //NEEDS TO BE NOT HARDCODED REL -- compare to that in the AlongStepDoIt
    //if( theDir.dot(surfaceNorm) <= epsilonDotProductForNorm ) continue;
    //We use the small tolerance in case there are floating point errors. It just needs to be much smaller than the deltaPhi
    if( theDir.dot(surfaceNorm) < (dotProductThreshold_Norm - smallTolerance) ) continue;
    //case we're talking about things that are being checked for being 86 degrees apart

    
    //In the case that the tangent vectors exist (i.e. we're in a "check for stuck QPs mode"), do a sanity check first
    if( tangVect1.mag() > 0 && tangVect2.mag() > 0 ){
      G4double minDot = tangVect1.dot(tangVect2); //Can move this outside loop to speed up REL?
      if( theDir.dot(tangVect1) < minDot || theDir.dot(tangVect2) < minDot ) continue;
    }

    //Now check safety
    G4double distToBound = 0;
    if( volIsMother ){ distToBound = theVolSolid->DistanceToOut(pos,theDir); }
    else{ distToBound = theVolSolid->DistanceToIn(pos,theDir); }
    if( distToBound < the2DSafety && distToBound > 0 ){ the2DSafety = distToBound; }    

    //Debugging
    if( verboseLevel > 5 ){
      G4cout << "C2DSFAB Function Point A | At angle: " << deltaPhi <<", (direction: " << theDir << "), distToBound: " << distToBound << ", dot: " << theDir.dot(surfaceNorm) << G4endl;
    }
    
  }

  //Rotate in the negative angular direction
  for( G4int iV = 0; iV < nVHalfCircle; ++iV ){
    G4double deltaPhi = -1.0*CLHEP::pi*((double)iV/(double)nVHalfCircle);
    theDir = surfaceNorm;
    theDir.rotateZ(deltaPhi); //Needs to be plane-agnostic at some point REL

    //G4double epsilonDotProductForNorm = 0.0896393089; //0.07; //NEEDS TO BE NOT HARDCODED REL -- compare to that in the AlongStepDoIt
    //if( theDir.dot(surfaceNorm) <= epsilonDotProductForNorm ) continue;
    if( theDir.dot(surfaceNorm) < (dotProductThreshold_Norm - smallTolerance) ) continue;
    
    //In the case that the tangent vectors exist (i.e. we're in a "check for stuck QPs mode"), do a sanity check first
    if( tangVect1.mag() > 0 && tangVect2.mag() > 0 ){
      G4double minDot = tangVect1.dot(tangVect2); //Can move this outside loop to speed up REL?
      if( theDir.dot(tangVect1) < minDot || theDir.dot(tangVect2) < minDot ) continue;
    }

    //Now check safety
    G4double distToBound = 0;
    if( volIsMother ){ distToBound = theVolSolid->DistanceToOut(pos,theDir); }
    else{ distToBound = theVolSolid->DistanceToIn(pos,theDir); }
    if( distToBound < the2DSafety && distToBound > 0 ){ the2DSafety = distToBound; }    

    //Debugging
    if( verboseLevel > 5 ){
      G4cout << "C2DSFAB Function Point B | At angle: " << deltaPhi <<", (direction: " << theDir << "), distToBound: " << distToBound << ", dot: " << theDir.dot(surfaceNorm) << G4endl;
    }
  }

  
  timestampEnd = clock();

  //Debugging
  if( verboseLevel > 5 ){
    G4cout << "C2DSFAB Function Point C | Time elapsed during 2D DistToOutLoop: " << double(timestampEnd-timestampStart)/double(CLOCKS_PER_SEC) << " seconds" << G4endl;
    G4cout << "C2DSFAB Function Point C | clocks_per_sec: " << CLOCKS_PER_SEC << G4endl;
    G4cout << "C2DSFAB Function Point C | In Compute2DSafetyFromABoundary, looking at a safety of: " << the2DSafety << ", and volIsMother is: " << volIsMother << G4endl;
  }
  return the2DSafety;
}

// Get normal to enclosing volume at boundary point in global coordinates
G4ThreeVector G4CMP::GetSurfaceNormal(const G4Step& step) {

  G4int verboseLevel = G4CMPConfigManager::GetVerboseLevel();
  if( verboseLevel > 5 ){
    G4cout << "---------- G4CMPGeometryUtils::GetSurfaceNormal ----------" << G4endl;
  }
  
  //RL: when we have nested geometries and the step is impinging upon an internal/daughter volume
  //from a parent volume, it seems that the returned normal is kind of wonky if one uses just the pre-step volume
  //as a way to access that. So we need a bit of logic to handle these situations.
  
  //First, let's get our pre-point volume touchable and its associated solid, and the same for the post-step
  const G4VTouchable* preTouch = step.GetPreStepPoint()->GetTouchable();
  const G4VTouchable* postTouch = step.GetPostStepPoint()->GetTouchable();
  G4VSolid* preSolid = preTouch->GetVolume()->GetLogicalVolume()->GetSolid();
  G4VSolid* postSolid = postTouch->GetVolume()->GetLogicalVolume()->GetSolid();

  //Check to see if the position is within spitting distance of either object's boundary. First, understand whether we're
  //inside either of the
  G4ThreeVector pos_prePV = step.GetPostStepPoint()->GetPosition();
  G4ThreeVector pos_postPV = step.GetPostStepPoint()->GetPosition();

  //Debugging
  if( verboseLevel > 5 ){
    G4cout << "GSN Function Point A | pos_prePV1: " << pos_prePV << ", in volume " << preTouch->GetVolume()->GetName() << G4endl;
    G4cout << "GSN Function Point A | pos_postPV1: " << pos_postPV << ", in volume " << postTouch->GetVolume()->GetName() << G4endl;
  }
    
  RotateToLocalPosition(preTouch, pos_prePV);
  RotateToLocalPosition(postTouch,pos_postPV);
  EInside postStepInPrePV = preSolid->Inside(pos_prePV);
  EInside postStepInPostPV = postSolid->Inside(pos_postPV);


  //Debugging
  if( verboseLevel > 5 ){
    G4cout << "GSN Function Point B | pos_prePV2: " << pos_prePV << G4endl;
    G4cout << "GSN Function Point B | pos_postPV2: " << pos_postPV << G4endl;  
    G4cout << "GSN Function Point B | PostStepInPrePV: " << postStepInPrePV << ", postStepInPostPV: " << postStepInPostPV << G4endl;
    G4cout << "GSN Function Point B | fabs(preSolid->DistanceToOut(pos_prePV)):" << fabs(preSolid->DistanceToOut(pos_prePV)) << ", fabs(preSolid->DistanceToIn(pos_prePV)): " << fabs(preSolid->DistanceToIn(pos_prePV)) << G4endl;
    G4cout << "GSN Function Point B | fabs(postSolid->DistanceToOut(pos_postPV)): " << fabs(postSolid->DistanceToOut(pos_postPV)) << ", fabs(postSolid->DistanceToIn(pos_postPV)): " << fabs(postSolid->DistanceToIn(pos_postPV)) << G4endl;
  }

  
  //Now, we have some logic. We don't want as much logic as in G4CMPBoundaryUtils::CheckBoundarySurface(), especially
  //when it comes to what happens when we don't have a point on a bonafide surface -- that code should be run
  //first to confirm that we're indeed on a boundary surface. But we do need something to tell us which surface's normal to reflect over.
  double tolerance = 1.0e-11 * mm; //Need to not hardcode this -- find a better way to implement. REL
  if( postStepInPrePV == kSurface ){
    //Reflect over the pre-PV normal
    G4ThreeVector preSolidNorm = preSolid->SurfaceNormal(pos_prePV);

    //Debugging
    if( verboseLevel > 5 ){ G4cout << "GSN Function Point C | returning pre-PV surface norm: " << preSolidNorm << G4endl; }
    RotateToGlobalDirection(preTouch,preSolidNorm);
    if( verboseLevel > 5 ){ G4cout << "GSN Function Point C | after rotation, this surface norm is: " << preSolidNorm << G4endl; }
    return preSolidNorm;
  }
  else if( postStepInPostPV == kSurface ){
    //Reflect over the post-PV normal
    G4ThreeVector postSolidNorm = postSolid->SurfaceNormal(pos_postPV);
    if( verboseLevel > 5 ){ G4cout << "GSN Function Point D | returning post-PV surface norm: " << postSolidNorm << G4endl; }
    RotateToGlobalDirection(postTouch,postSolidNorm);
    if( verboseLevel > 5 ){ G4cout << "GSN Function Point D | after rotation, this surface norm is: " << postSolidNorm << G4endl; }
    return postSolidNorm;
  }
  //If we are very near a surface and within tolerance, still okay -- this is a pre-PV reflection
  else if( (fabs(preSolid->DistanceToOut(pos_prePV)) > 0 && fabs(preSolid->DistanceToOut(pos_prePV)) < tolerance ) ||
	   (fabs(preSolid->DistanceToIn(pos_prePV)) > 0 && fabs(preSolid->DistanceToIn(pos_prePV)) < tolerance ) ){
    G4ThreeVector preSolidNorm = preSolid->SurfaceNormal(pos_prePV);
    if( verboseLevel > 5 ){ G4cout << "GSN Function Point E | returning pre-PV surface norm: " << preSolidNorm << G4endl; }
    RotateToGlobalDirection(preTouch,preSolidNorm);
    if( verboseLevel > 5 ){ G4cout << "GSN Function Point E | after rotation, this surface norm is: " << preSolidNorm << G4endl; }
    return preSolidNorm;
  }
  //If we are very near a surface and within tolerance, still okay -- this is for post-PV reflection
  else if( (fabs(postSolid->DistanceToOut(pos_postPV)) > 0 && fabs(postSolid->DistanceToOut(pos_postPV)) < tolerance ) ||
	   (fabs(postSolid->DistanceToIn(pos_postPV)) > 0 && fabs(postSolid->DistanceToIn(pos_postPV)) < tolerance ) ){
    G4ThreeVector postSolidNorm = postSolid->SurfaceNormal(pos_postPV);
    if( verboseLevel > 5 ){ G4cout << "GSN Function Point F | returning post-PV surface norm: " << postSolidNorm << G4endl; }
    RotateToGlobalDirection(postTouch,postSolidNorm);
    if( verboseLevel > 5 ){ G4cout << "GSN Function Point F | after rotation, this surface norm is: " << postSolidNorm << G4endl; }
    return postSolidNorm;
  }
  //Otherwise, we're not on (or within tolerance of) a surface.
  else{
    
    //Throw an error -- shouldn't ever be here since it means we're not on a surface
    G4ExceptionDescription msg;
    msg << "G4CMP::GetSurfaceNormal() seems to think we're not on a surface. Volumes are: "
	<< preSolid->GetName() << " and " << postSolid->GetName() << G4endl;
    G4Exception("G4CMP::GetSurfaceNormal()", "Geometry00X",
		FatalException, msg);
    G4ThreeVector dummy(0,0,0);
    return dummy;    
  }
}



// Create non-tracking Navigator for use with position finding below
G4Navigator* G4CMP::GetNavigator() {
  static G4ThreadLocal G4Navigator* theNavigator = 0;
  if (!theNavigator) theNavigator = new G4Navigator;

  // Make sure current world volume is the one in use  
  G4VPhysicalVolume* theWorld =
    G4TransportationManager::GetTransportationManager()->
      GetNavigatorForTracking()->GetWorldVolume();

  if (theNavigator->GetWorldVolume() != theWorld)
    theNavigator->SetWorldVolume(theWorld);

  return theNavigator;
}

// Get placement volume at specified global position

G4VPhysicalVolume* G4CMP::GetVolumeAtPoint(const G4ThreeVector& pos) {
  return GetNavigator()->LocateGlobalPointAndSetup(pos,0,false);
}


// Get touchable (geometry tree path) for specified global position

G4VTouchable* G4CMP::CreateTouchableAtPoint(const G4ThreeVector& pos) {
  G4VTouchable* touchable = new G4TouchableHistory;

  GetNavigator()->LocateGlobalPointAndUpdateTouchable(pos, touchable, false);

  // Sanity check: touchable's volume should match GetVolumeAtPoint()
#ifdef G4CMP_DEBUG
  if (touchable->GetVolume() != GetVolumeAtPoint(pos)) {
    G4ExceptionDescription msg;
    msg << "Position " << pos << " returns different results for touchable"
	<< " (" << touchable->GetVolume()->GetName() << ") vs. GetVolume"
	<< " (" << GetVolumeAtPoint(pos)->GetName() << ")" << G4endl;
    G4Exception("G4CMP::CreateTouchableAtPoint", "Geometry008",
		EventMustBeAborted, msg);
  }
#endif

  return touchable;
}


// Adjust specified position to avoid surface of current (touchable) volume

G4ThreeVector G4CMP::ApplySurfaceClearance(const G4VTouchable* touch,
					   G4ThreeVector pos) {
  // Clearance is the minimum distance where a position is guaranteed Inside
  const G4double clearance = G4CMPConfigManager::GetSurfaceClearance();

  // If the step is near a boundary, create the secondary in the initial volume
  G4VPhysicalVolume* pv = touch->GetVolume();
  G4ThreadLocalStatic auto latMan = G4LatticeManager::GetLatticeManager();
  G4LatticePhysical* lat = latMan->GetLattice(pv);

  if (!lat) {		// No lattice in touchable's volume, try pos instead
    pv = G4CMP::GetVolumeAtPoint(pos);
    lat = latMan->GetLattice(pv);

    if (!lat) {
      G4ExceptionDescription msg;
      msg << "Position " << pos << " not associated with valid volume.";
      G4Exception("G4CMP::CreateSecondary", "Secondary008",
		  EventMustBeAborted, msg);
      return pos;
    }
  }

  // Work in local coordinates, adjusting position to be clear of surface
  RotateToLocalPosition(touch, pos);

  G4VSolid* solid = pv->GetLogicalVolume()->GetSolid();
  G4ThreeVector norm = solid->SurfaceNormal(pos);

  while (solid->Inside(pos) != kInside ||
	 solid->DistanceToOut(pos,norm) < clearance) {
    pos -= norm*clearance;
    norm = solid->SurfaceNormal(pos);	// Nearest surface may change
  }

  RotateToGlobalPosition(touch, pos);
  return pos;
}

G4double G4CMP::ComputeDotProductThreshold_Norm(int full_circle_nV)
{
  G4int half_circle_nV = full_circle_nV / 2;
  if( half_circle_nV % 2 != 0 ){
    G4ExceptionDescription msg;
    msg << "G4CMP::ComputeDotProductThreshold_Norm seems to see that half_circle_nV is not divisible by two. Please fix." << G4endl;
    G4Exception("G4CMP::ComputeDotProductThreshold_Norm()", "Geometry00X",FatalException, msg);
  }

  
  //Divide the nV into 180 to understand the degrees per step
  G4double degreesPerStep = 180.0/((double)half_circle_nV);

  
  //Since half_circle_nV should be divisible by 2, we should always have a step at exactly 90 degrees.
  //We want to identify the angle associated with one step prior to 90 degrees
  G4double preTangentStep_deg = 90.0-degreesPerStep;  
  G4double dotProductThreshold_norm = cos(CLHEP::pi/180.*preTangentStep_deg);
  return dotProductThreshold_norm;
}


G4double G4CMP::ComputeDotProductThreshold_Tang(int full_circle_nV)
{
  G4int half_circle_nV = full_circle_nV / 2;
  if( half_circle_nV % 2 != 0 ){
    G4ExceptionDescription msg;
    msg << "G4CMP::ComputeDotProductThreshold_Tang seems to see that half_circle_nV is not divisible by two. Please fix." << G4endl;
    G4Exception("G4CMP::ComputeDotProductThreshold_Tang()", "Geometry00X",FatalException, msg);
  }
  
  //Divide the nV into 180 to understand the degrees per step
  G4double degreesPerStep = 180.0/((double)half_circle_nV);

  //Since half_circle_nV should be divisible by 2, we should always have a step at exactly 90 degrees.
  //We want to identify the angle associated with one step prior to 90 degrees
  G4double dotProductThreshold_tang = cos(CLHEP::pi/180.*degreesPerStep);
  return dotProductThreshold_tang;
}

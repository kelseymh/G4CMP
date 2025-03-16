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

/*
//Find safety in XY, if we're ON a boundary. Here, since we're on a boundary, the DistToIn and DistToOut
//functions will be pretty useless (will return 0), but since hopefully boundary interactions won't
//actually occur that frequently, we can just do the spamming of vectors and find the smallest nonzero
//value. We'll also want to watch out for curvature, in a scenario where the random direction is quite
//parallel to the boundary. Will need to debug this later.
G4double G4CMP::Get2DSafetyFromBoundary(const G4Track& theTrack)
{
  G4cout << "REL in Get2DSafetyFromBoundary" << G4endl;
  
  //First, identify the touchable. Recall that for use in the GPIL functions, we don't have a post-step point yet, so the
  //touchable should correspond to the pre-step point. Identify the volume name, get 
  const G4VTouchable* volTouch = theTrack.GetStep()->GetPreStepPoint()->GetTouchable();
  G4VPhysicalVolume* volPhys = volTouch->GetVolume();
  G4LogicalVolume* volLog = volPhys->GetLogicalVolume();
  G4VSolid * volSolid = volLog->GetSolid();
  G4ThreeVector pos = theTrack.GetPosition();

  //Rotate the position to the touchable's coordinates
  RotateToLocalPosition(volTouch, pos);

  //Spam a set of DistToOuts in different directions. This part slows things down substantially and should be used sparingly.
  //This part should be sped up a bit by figuring out the normal and then only throwing vectors that are on the "good" side
  //of that normal
  double nearestSurfaceDistXY = DBL_MAX;
  G4ThreeVector theDir;
  G4int nV = 100;
  clock_t timestampStart, timestampEnd;
  timestampStart = clock();
  for( G4int iV = 0; iV < nV; ++iV ){
    G4double phi = 2*CLHEP::pi*((double)iV/(double)nV);
    G4double x = cos(phi);
    G4double y = sin(phi);
    theDir.setX(x);
    theDir.setY(y);
    theDir.setZ(0);
    G4double distToOut = volSolid->DistanceToOut(pos,theDir);
    G4cout << "---> REL/EY: Running hack-y 2-D safety. For phi of " << phi*180/CLHEP::pi << " deg, DistToOut returns a dist to boundary of: " << distToOut << G4endl;
    if( distToOut < nearestSurfaceDistXY && distToOut > 0 ) nearestSurfaceDistXY = distToOut;
  }
  timestampEnd = clock();
  G4cout << "Time elapsed during 2D DistToOutLoop: " << double(timestampEnd-timestampStart)/double(CLOCKS_PER_SEC) << " seconds" << G4endl;
  G4cout << "clocks_per_sec: " << CLOCKS_PER_SEC << G4endl;
  
  
  //Loop through the daughters of this mother volume and find distances to them
  for( int iD = 0; iD < volLog->GetNoDaughters(); ++iD ){
    G4VPhysicalVolume * volDaughterPhys = volLog->GetDaughter(iD);
    
    //This stuff is ripped verbatim from G4NormalNavigation. We're going to see how well it works...
    G4AffineTransform sampleTf(volDaughterPhys->GetRotation(),
			       volDaughterPhys->GetTranslation());

    sampleTf.Invert();
    const G4ThreeVector samplePoint = sampleTf.TransformPoint(pos);
    const G4VSolid * volDaughterSolid = volDaughterPhys->GetLogicalVolume()->GetSolid();
    const G4double volDaughterSafety = volDaughterSolid->DistanceToIn(samplePoint);

    //Since we're *outside* the daughter volume laterally, it means that the distToIn is in fact the distance in XY. No additional
    //math is needed.
    if( volDaughterSafety < nearestSurfaceDistXY ) nearestSurfaceDistXY = volDaughterSafety;
    
    //For debugging, cout
    G4cout << "REL: In Get2DSafetyFromBoundary, we are looking at a daughter volume of: " << volDaughterPhys->GetName() << " which has a distToIn of " << volDaughterSafety << G4endl;
    
  }

  G4cout << "REL: In Get2DSafetyFromBoundary, returning nearestSurfaceDistXY: " << nearestSurfaceDistXY << G4endl;
  return nearestSurfaceDistXY; 
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
    G4cout << "Overall safety during mother safety check: " << overallSafety << G4endl;
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
    G4cout << "Overall safety during daughter safety check " << iD << ": " << overallSafety << G4endl;
  }
  
  //Return the safety
  return overallSafety; 
  
}

//Looking "inward" within a mother volume to its daughter volumes to identify safeties. Generally, don't want to use this on its
//own. Should only be called from Get2DSafety, not by separate classes.
G4double G4CMP::Compute2DSafetyToDaughterVolume(const G4ThreeVector & pos, const G4ThreeVector & momDir, G4LogicalVolume * motherLog, bool safetyFromABoundary, G4int daughterID, G4ThreeVector surfaceNorm, G4ThreeVector tangVect1, G4ThreeVector tangVect2 ){

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
  

  /* Old, simpler logic
  if( safetyFromABoundary ){
    if( volDaughterSafety < safety && volDaughterSafety != 0 ) safety = volDaughterSafety;
  }
  else{
    if( volDaughterSafety < safety ) safety = volDaughterSafety;
  }
  */

  //More complicated logic, which takes into account the effect of being on a boundary of a surface with
  //concavities.
  if( safetyFromABoundary ){  
    G4cout << "SafetyFromABoundary is triggering, with volDaughterSafety: " << volDaughterSafety << G4endl;

    //Determine if the current boundary belongs to this daughter
    G4bool currentBoundaryBelongsToThisDaughter = false;
    G4double boundaryTolerance = 1e-12; //REL HARDCODED NEED TO FIX -- maybe use kSurface here instead?
    if( volDaughterSafety < boundaryTolerance ){
      currentBoundaryBelongsToThisDaughter = true;
      G4cout << "Changing currentBoundaryBelongsToThisDaughter to true." << G4endl;
    }
    
    //If it doesn't, then we can just proceed as normal -- I think there aren't odd edge cases here
    if( !currentBoundaryBelongsToThisDaughter ){
      G4cout << "We're running the standard thing " << G4endl;
      if( volDaughterSafety < safety ) safety = volDaughterSafety;
    }    
    //If it does, we have do do a sweep so we can see if there are places on this daughter volume
    //that are closer to the given point than the safeties to other daughter volumes (and the mother)
    else{
      G4cout << "Somehow we think that this current boundary belongs to this daughter." << G4endl;
      G4double safetyToThisDaughterBoundary = Compute2DSafetyToThisDaughterBoundary(volDaughterSolid,samplePoint,rotatedSurfaceNorm,rotatedTangVect1,rotatedTangVect2);
      if( safetyToThisDaughterBoundary < safety ) safety = safetyToThisDaughterBoundary;
    }
  }
  //If we're not on a boundary we get easy logic
  else{
    if( volDaughterSafety < safety ) safety = volDaughterSafety;
  }
  
    
  
  /* //THIS WILL NEED TO GO IN SOMEHOW BUT ONCE WE CAN CONFIRM THINGS ARE WORKING WITH SIMPLER GEOMETRY FIRST
  //Handle the logic for whether we're trying to compute from a boundary. If we're on a boundary, we'll want to check
  //if:
  //1. That boundary belongs to this daughter (safety = 0)
  //   a. This daughter also possesses the next-closest boundary
  //   b. This daughter does not possess the next closest boundary
  //2. That boundary does not belong to this daughter (safety != 0 )
  //It's hard to differentiate points 1a and 1b without doing a search in direction space.
  //So we'll do that if the distToIn here is zero.
  if( safetyFromABoundary ){

    //We're ON this daughter boundary. Now do a scan.
    if( volDaughterSafety == 0 ){

      //First, get a surface normal at this point. 
      G4ThreeVector surfaceNorm = volDaughterSolid->SurfaceNormal(samplePoint);      
      
      //Second, use the surface normal with the momentum direction (which should be in-plane with the film)
      //to define a plane in which we will do a directional scan.
      //This actually bypasses the need to hardcode which plane/dimensions the scan is in. First, check the direction
      /
      
      
      
      
    }

    
    //At this point, there are two scenarios we should be aware of:
    //1. We're on a mother boundary, and all distances to daughters are nonzero.
    //2. We're on a daughter boundary. In this case, there are two sub-cases
    //   a. We're on a daughter boundary and the next-closest boundary in the half-circle of the surface normal
    //      is on a different daughter 
    //   b. We're on a daughter boundary and the next-closest boundary in the half-circle of the surface normal
    //      is on the same daughter.
    //It's hard to identify whether we're in scenario 2a or 2b without doing a sweep of direction vectors. So to be
    //safe, if this daughter has a "simple" safety of 0 (i.e. we're on its surface), we're going to compute the
    //DistanceToIn for itself over a range of directions. This range of directions cannot span a perfect 180 degrees
    //along the surface because the space may be convex, so we'll have it span 180-2*eta, where eta is a small angle
    //into which we'll never end up sending the particles.

    //First, check to see if volDaughterSafety is less than the 
    if( volDaughterSafety < safety )


    //The logic for 1 and 2 are the same: compute a simple distance to all daughters, and ignore any for which
    //that distanceToIn is zero (that's the one we're on).
    //The logic for 3 is significantly more complicated, especially in the scenario where we have a convex space. So we'll
    //want to do something to help with that. 
    
    if( volDaughterSafety < safety && volDaughterSafety > 0 ) safety = volDaughterSafety;
  }
  else{
    if( volDaughterSafety < safety ) safety = volDaughterSafety;
  }
  */



  G4cout << "REL: In Get2DSafetyToDaughterVolume, we are looking at a daughter volume of: " << volDaughterPhys->GetName() << " which has a distToIn of " << safety << G4endl;    
  return safety;
}

//This is a specific 2D safety computation, from a daughter boundary to its own boundary. Only really useful if the surface has a concavity
//that we're looking into. This has code overlapping with the mother version of this but for generality in the long-term, it may be more
//organized to split into "look for 2D safety from daughters" and "look for 2D safty to mother"
G4double G4CMP::Compute2DSafetyToThisDaughterBoundary(const G4VSolid * volDaughterSolid,
						      G4ThreeVector samplePoint,
						      G4ThreeVector rotatedSurfaceNorm,
						      G4ThreeVector rotatedTangVect1,
						      G4ThreeVector rotatedTangVect2)
{
  //Debugging
  G4cout << "---------- G4CMPGeometryUtils::Compute2DSafetyToThisDaughterBoundary() ----------" << G4endl;
  G4cout << "C2DSTTB Function Point A | Rotated surface norm: " << rotatedSurfaceNorm << G4endl;

  //Spam a set of DistToOuts in different directions. This part slows things down substantially and should be used sparingly. This part
  //should be replaced with geometry math
  double safety = DBL_MAX;
  G4ThreeVector theDir;
  clock_t timestampStart, timestampEnd;
  timestampStart = clock();

  //Want to generate our points such that we start from the norm and work outwards
  G4int nV = 70; //Parameterize this as half of the one used for the bulk -- we'll sweep both directions
  
  //Rotate in positive angular direction
  for( G4int iV = 0; iV < nV; ++iV ){
    G4double deltaPhi = 1.0*CLHEP::pi*((double)iV/(double)nV);
    theDir = rotatedSurfaceNorm;
    theDir.rotateZ(deltaPhi); //Needs to be plane-agnostic at some point REL

    G4double epsilonDotProductForNorm = 0.0896393089; //Used to be 0.07//NEEDS TO BE NOT HARDCODED REL -- compare to that in the AlongStepDoIt
    if( theDir.dot(rotatedSurfaceNorm) <= epsilonDotProductForNorm ) continue;
    
    //In the case that the tangent vectors exist (i.e. we're in a "check for stuck QPs mode"), do a sanity check first
    if( rotatedTangVect1.mag() > 0 && rotatedTangVect2.mag() > 0 ){
      G4double minDot = rotatedTangVect1.dot(rotatedTangVect2); //Can move this outside loop to speed up REL?
      if( theDir.dot(rotatedTangVect1) < minDot || theDir.dot(rotatedTangVect2) < minDot ) continue;
    }

    //Now check safety
    G4double thisDaughterDirectionalDistToIn = volDaughterSolid->DistanceToIn(samplePoint,theDir);
    if( thisDaughterDirectionalDistToIn < safety && thisDaughterDirectionalDistToIn > 0 ) safety = thisDaughterDirectionalDistToIn;

    //Debugging
    G4cout << "C2DSTTDB Function Point B | At angle: " << deltaPhi << ", (direction: " << theDir << "), thisDaughterDirectionalDistToIn: " << thisDaughterDirectionalDistToIn << G4endl;
  }

  //Rotate in the negative angular direction
  for( G4int iV = 0; iV < nV; ++iV ){
    G4double deltaPhi = -1.0*CLHEP::pi*((double)iV/(double)nV);
    theDir = rotatedSurfaceNorm;
    theDir.rotateZ(deltaPhi); //Needs to be plane-agnostic at some point REL

    G4double epsilonDotProductForNorm = 0.0896393089; //0.07; //NEEDS TO BE NOT HARDCODED REL -- compare to that in the AlongStepDoIt
    if( theDir.dot(rotatedSurfaceNorm) <= epsilonDotProductForNorm ) continue;
    
    //In the case that the tangent vectors exist (i.e. we're in a "check for stuck QPs mode"), do a sanity check first
    if( rotatedTangVect1.mag() > 0 && rotatedTangVect2.mag() > 0 ){
      G4double minDot = rotatedTangVect1.dot(rotatedTangVect2); //Can move this outside loop to speed up REL?
      if( theDir.dot(rotatedTangVect1) < minDot || theDir.dot(rotatedTangVect2) < minDot ) continue;
    }

    //Now check safety
    G4double thisDaughterDirectionalDistToIn = volDaughterSolid->DistanceToIn(samplePoint,theDir);
    if( thisDaughterDirectionalDistToIn < safety && thisDaughterDirectionalDistToIn > 0 ) safety = thisDaughterDirectionalDistToIn;

    //Debugging
    G4cout << "C2DSTTDB Function Point C | At angle: " << deltaPhi <<", (direction: " << theDir << "), thisDaughterDirectionalDistToIn: " << thisDaughterDirectionalDistToIn << G4endl;
  }
  timestampEnd = clock();
  G4cout << "In Compute2DSafetyToThisDaughterBoundary, looking at a safety of: " << safety << G4endl;
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
    G4double motherSafetyFromABoundary = Compute2DMotherSafetyFromABoundary(motherSolid,pos,surfaceNorm,tangVect1,tangVect2);
    motherSafety = motherSafetyFromABoundary;
  }
  //Otherwise, doesn't matter -- this is freeform
  else{
    G4double motherSafetyFromTheBulk = Compute2DMotherSafetyFromtheBulk(motherSolid,pos);
    motherSafety = motherSafetyFromTheBulk;
  }


  /*
  //Last, there is an edge case that may arise if we're not fine enough with our sampling, which arises if we try to
  //take a step basically right along the tangent vectors. (We note that this "spam in several directions" strategy may
  //break down if a surface is not "simple" in phi, but this is a scenario where even "simple" geometries break down.) For
  //this, we should re-compute the safety along exactly the tangent vectors
  if( tangVect1.mag() > 0 && tangVect2.mag() > 0 ){
    G4double safetyTang1 = motherSolid->DistanceToOut(pos,tangVect1);
    G4double safetyTang2 = motherSolid->DistanceToOut(pos,tangVect2);
    G4cout << "Constrained Safety along Tang1: " << safetyTang1 << ", constrained safety along Tang2: " << safetyTang2 << G4endl;
    if( safetyTang1 < motherSafety ) motherSafety = safetyTang1;
    if( safetyTang2 < motherSafety ) motherSafety = safetyTang2;
  }
  */
  return motherSafety;
}

//Compute the 2D safety to the mother from a point in the bulk (i.e. not on a boundary)
G4double G4CMP::Compute2DMotherSafetyFromtheBulk(const G4VSolid * motherSolid,G4ThreeVector pos)
{
  //Spam a set of DistToOuts in different directions. This part slows things down substantially and should be used sparingly. This part
  //should be replaced with geometry math
  double motherSafety = DBL_MAX;
  G4ThreeVector theDir;
  G4int nV = 140;
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
  G4cout << "Time elapsed during 2D DistToOutLoop: " << double(timestampEnd-timestampStart)/double(CLOCKS_PER_SEC) << " seconds" << G4endl;
  G4cout << "clocks_per_sec: " << CLOCKS_PER_SEC << G4endl;
  G4cout << "In Compute2DMotherSafetyFromTheBulk, looking at a mother safety of: " << motherSafety << G4endl;
  return motherSafety;
}

//Compute the 2D safety to the mother from a point on the surface. The vectors here are all in the mother
//frame, and the scan happens in both directions to cleanly define a set of absolute angles with respect to the tangents
G4double G4CMP::Compute2DMotherSafetyFromABoundary(const G4VSolid * motherSolid,G4ThreeVector pos, G4ThreeVector surfaceNorm, G4ThreeVector tangVect1, G4ThreeVector tangVect2)
{
  //Spam a set of DistToOuts in different directions. This part slows things down substantially and should be used sparingly. This part
  //should be replaced with geometry math
  double motherSafety = DBL_MAX;
  G4ThreeVector theDir;
  clock_t timestampStart, timestampEnd;
  timestampStart = clock();

  //Want to generate our points such that we start from the norm and work outwards
  G4int nV = 70; //Parameterize this as half of the one used for the bulk -- we'll sweep both directions
  
  //Rotate in positive angular direction
  for( G4int iV = 0; iV < nV; ++iV ){
    G4double deltaPhi = 1.0*CLHEP::pi*((double)iV/(double)nV);
    theDir = surfaceNorm;
    theDir.rotateZ(deltaPhi); //Needs to be plane-agnostic at some point REL

    //Check to make sure that the scan vector is sufficiently away from the surface tangent so we don't end up with
    //ridiculously crazy small steps'
    
    G4double epsilonDotProductForNorm = 0.0896393089; //0.07; //NEEDS TO BE NOT HARDCODED REL -- compare to that in the AlongStepDoIt
    if( theDir.dot(surfaceNorm) <= epsilonDotProductForNorm ) continue;
    //case we're talking about things that are being checked for being 86 degrees apart

    
    //In the case that the tangent vectors exist (i.e. we're in a "check for stuck QPs mode"), do a sanity check first
    if( tangVect1.mag() > 0 && tangVect2.mag() > 0 ){
      G4double minDot = tangVect1.dot(tangVect2); //Can move this outside loop to speed up REL?
      if( theDir.dot(tangVect1) < minDot || theDir.dot(tangVect2) < minDot ) continue;
    }

    //Now check safety
    G4double distToOut = motherSolid->DistanceToOut(pos,theDir);
    if( distToOut < motherSafety && distToOut > 0 ){ motherSafety = distToOut; }

    //Debugging
    G4cout << "C2DMSFAB Function Point A | At angle: " << deltaPhi <<", (direction: " << theDir << "), distToOut: " << distToOut << ", dot: " << theDir.dot(surfaceNorm) << G4endl;
    
  }

  //Rotate in the negative angular direction
  for( G4int iV = 0; iV < nV; ++iV ){
    G4double deltaPhi = -1.0*CLHEP::pi*((double)iV/(double)nV);
    theDir = surfaceNorm;
    theDir.rotateZ(deltaPhi); //Needs to be plane-agnostic at some point REL

    G4double epsilonDotProductForNorm = 0.0896393089; //0.07; //NEEDS TO BE NOT HARDCODED REL -- compare to that in the AlongStepDoIt
    if( theDir.dot(surfaceNorm) <= epsilonDotProductForNorm ) continue;
    
    //In the case that the tangent vectors exist (i.e. we're in a "check for stuck QPs mode"), do a sanity check first
    if( tangVect1.mag() > 0 && tangVect2.mag() > 0 ){
      G4double minDot = tangVect1.dot(tangVect2); //Can move this outside loop to speed up REL?
      if( theDir.dot(tangVect1) < minDot || theDir.dot(tangVect2) < minDot ) continue;
    }

    //Now check safety
    G4double distToOut = motherSolid->DistanceToOut(pos,theDir);
    if( distToOut < motherSafety && distToOut > 0 ){ motherSafety = distToOut; }

    //Debugging
    G4cout << "C2DMSFAB Function Point B | At angle: " << deltaPhi <<", (direction: " << theDir << "), distToOut: " << distToOut << ", dot: " << theDir.dot(surfaceNorm) << G4endl;
  }

  timestampEnd = clock();
  G4cout << "Time elapsed during 2D DistToOutLoop: " << double(timestampEnd-timestampStart)/double(CLOCKS_PER_SEC) << " seconds" << G4endl;
  G4cout << "clocks_per_sec: " << CLOCKS_PER_SEC << G4endl;
  G4cout << "In Compute2DMotherSafetyFromABoundary, looking at a mother safety of: " << motherSafety << G4endl;
  return motherSafety;
}

/*  
//Hacky way to try to find safety in XY
G4double G4CMP::Get2DSafety(const G4Track& theTrack, G4ThreeVector & directionToNearestBoundary )
{
  G4cout << "REL in Get2DSafety" << G4endl;
  
  //The "final" strategy here will be a hybrid approach. For "easy" geometrical objects that are likely to be
  //used in thin film geometries, we'll do analytical calculations which will speed up the code execution. For now
  //we are going to do the slow thing just to get a first pass-through of the full physics. We'll come back and improve
  //speed in a bit.
  
  //First, identify the touchable. Recall that for use in the GPIL functions, we don't have a post-step point yet, so the
  //touchable should correspond to the pre-step point. Identify the volume name, get 
  const G4VTouchable* volTouch = theTrack.GetStep()->GetPreStepPoint()->GetTouchable();
  G4VPhysicalVolume* volPhys = volTouch->GetVolume();
  G4LogicalVolume* volLog = volPhys->GetLogicalVolume();
  G4VSolid * volSolid = volLog->GetSolid();
  G4ThreeVector pos = theTrack.GetPosition();

  //Rotate the position to the touchable's coordinates
  RotateToLocalPosition(volTouch, pos);

  //Spam a set of DistToOuts in different directions. This part slows things down substantially and should be used sparingly. This part
  //should be replaced with geometry math
  double nearestSurfaceDistXY = DBL_MAX;
  directionToNearestBoundary.setX(DBL_MAX);
  directionToNearestBoundary.setY(DBL_MAX);
  directionToNearestBoundary.setZ(DBL_MAX);  
  G4ThreeVector theDir;
  G4int nV = 100;
  clock_t timestampStart, timestampEnd;
  timestampStart = clock();
  for( G4int iV = 0; iV < nV; ++iV ){
    G4double phi = 2*CLHEP::pi*((double)iV/(double)nV);
    G4double x = cos(phi);
    G4double y = sin(phi);
    theDir.setX(x);
    theDir.setY(y);
    theDir.setZ(0);
    G4double distToOut = volSolid->DistanceToOut(pos,theDir);
    G4cout << "---> REL/EY: Running hack-y 2-D safety. For phi of " << phi*180/CLHEP::pi << " deg, DistToOut returns a dist to boundary of: " << distToOut << G4endl;
    if( distToOut < nearestSurfaceDistXY ){
      nearestSurfaceDistXY = distToOut;
      directionToNearestBoundary = theDir;
    }
  }
  timestampEnd = clock();
  G4cout << "Time elapsed during 2D DistToOutLoop: " << double(timestampEnd-timestampStart)/double(CLOCKS_PER_SEC) << " seconds" << G4endl;
  G4cout << "clocks_per_sec: " << CLOCKS_PER_SEC << G4endl;
  
  
  //Loop through the daughters of this mother volume
  int closestDaughter = -1;
  for( int iD = 0; iD < volLog->GetNoDaughters(); ++iD ){
    G4VPhysicalVolume * volDaughterPhys = volLog->GetDaughter(iD);
    
    //This stuff is ripped verbatim from G4NormalNavigation. We're going to see how well it works...
    G4AffineTransform sampleTf(volDaughterPhys->GetRotation(),
			       volDaughterPhys->GetTranslation());

    sampleTf.Invert();
    const G4ThreeVector samplePoint = sampleTf.TransformPoint(pos);
    const G4VSolid * volDaughterSolid = volDaughterPhys->GetLogicalVolume()->GetSolid();
    const G4double volDaughterSafety = volDaughterSolid->DistanceToIn(samplePoint);

    //Since we're *outside* the daughter volume laterally, it means that the distToIn is in fact the distance in XY. No additional
    //math is needed.
    if( volDaughterSafety < nearestSurfaceDistXY ){
      nearestSurfaceDistXY = volDaughterSafety;
      closestDaughter = iD;
    }
    
    //For debugging, cout
    G4cout << "REL: In Get2DSafety, we are looking at a daughter volume of: " << volDaughterPhys->GetName() << " which has a distToIn of " << volDaughterSafety << G4endl;

  }

  
  

  G4cout << "REL: In Get2DSafety, returning nearestSurfaceDistXY: " << nearestSurfaceDistXY << G4endl;
  return nearestSurfaceDistXY; 
}
*/

// Get normal to enclosing volume at boundary point in global coordinates
G4ThreeVector G4CMP::GetSurfaceNormal(const G4Step& step) {

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
  G4cout << "pos_prePV1: " << pos_prePV << ", in volume " << preTouch->GetVolume()->GetName() << G4endl;
  G4cout << "pos_postPV1: " << pos_postPV << ", in volume " << postTouch->GetVolume()->GetName() << G4endl;
  
  RotateToLocalPosition(preTouch, pos_prePV);
  RotateToLocalPosition(postTouch,pos_postPV);
  EInside postStepInPrePV = preSolid->Inside(pos_prePV);
  EInside postStepInPostPV = postSolid->Inside(pos_postPV);
  G4cout << "pos_prePV2: " << pos_prePV << G4endl;
  G4cout << "pos_postPV2: " << pos_postPV << G4endl;

  
  G4cout << "PostStepInPrePV: " << postStepInPrePV << ", postStepInPostPV: " << postStepInPostPV << ", fabs(preSolid->DistanceToOut(pos_prePV)):" << fabs(preSolid->DistanceToOut(pos_prePV)) << ", fabs(preSolid->DistanceToIn(pos_prePV)): " << fabs(preSolid->DistanceToIn(pos_prePV)) << ", fabs(postSolid->DistanceToOut(pos_postPV)): " << fabs(postSolid->DistanceToOut(pos_postPV)) << ", fabs(postSolid->DistanceToIn(pos_postPV)): " << fabs(postSolid->DistanceToIn(pos_postPV)) << G4endl;

  
  //Now, we have some logic. We don't want as much logic as in G4CMPBoundaryUtils::CheckBoundarySurface(), especially
  //when it comes to what happens when we don't have a point on a bonafide surface -- that code should be run
  //first to confirm that we're indeed on a boundary surface. But we do need something to tell us which surface's normal to reflect over.
  double tolerance = 1.0e-11 * mm; //Need to not hardcode this -- find a better way to implement. REL
  if( postStepInPrePV == kSurface ){
    //Reflect over the pre-PV normal
    G4ThreeVector preSolidNorm = preSolid->SurfaceNormal(pos_prePV);
    G4cout << "REL in G4CMPGeometryUtils::GetSurfaceNormal() --- returning pre-PV surface norm: " << preSolidNorm << G4endl;
    RotateToGlobalDirection(preTouch,preSolidNorm);
    G4cout << "REL in G4CMPGeometryUtils::GetSurfaceNormal() --- after rotation, this surface norm is: " << preSolidNorm << G4endl;
    return preSolidNorm;
  }
  else if( postStepInPostPV == kSurface ){
    //Reflect over the post-PV normal
    G4ThreeVector postSolidNorm = postSolid->SurfaceNormal(pos_postPV);
    G4cout << "REL in G4CMPGeometryUtils::GetSurfaceNormal() --- returning post-PV surface norm: " << postSolidNorm << G4endl;
    RotateToGlobalDirection(postTouch,postSolidNorm);
    G4cout << "REL in G4CMPGeometryUtils::GetSurfaceNormal() --- after rotation, this surface norm is: " << postSolidNorm << G4endl;
    return postSolidNorm;
  }
  //If we are very near a surface and within tolerance, still okay -- this is a pre-PV reflection
  else if( (fabs(preSolid->DistanceToOut(pos_prePV)) > 0 && fabs(preSolid->DistanceToOut(pos_prePV)) < tolerance ) ||
	   (fabs(preSolid->DistanceToIn(pos_prePV)) > 0 && fabs(preSolid->DistanceToIn(pos_prePV)) < tolerance ) ){
    G4ThreeVector preSolidNorm = preSolid->SurfaceNormal(pos_prePV);
    G4cout << "REL in G4CMPGeometryUtils::GetSurfaceNormal() --- returning pre-PV surface norm: " << preSolidNorm << G4endl;
    RotateToGlobalDirection(preTouch,preSolidNorm);
    G4cout << "REL in G4CMPGeometryUtils::GetSurfaceNormal() --- after rotation, this surface norm is: " << preSolidNorm << G4endl;
    return preSolidNorm;
  }
  //If we are very near a surface and within tolerance, still okay -- this is for post-PV reflection
  else if( (fabs(postSolid->DistanceToOut(pos_postPV)) > 0 && fabs(postSolid->DistanceToOut(pos_postPV)) < tolerance ) ||
	   (fabs(postSolid->DistanceToIn(pos_postPV)) > 0 && fabs(postSolid->DistanceToIn(pos_postPV)) < tolerance ) ){
    G4ThreeVector postSolidNorm = postSolid->SurfaceNormal(pos_postPV);
    G4cout << "REL in G4CMPGeometryUtils::GetSurfaceNormal() --- returning post-PV surface norm: " << postSolidNorm << G4endl;
    RotateToGlobalDirection(postTouch,postSolidNorm);
    G4cout << "REL in G4CMPGeometryUtils::GetSurfaceNormal() --- after rotation, this surface norm is: " << postSolidNorm << G4endl;
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



/*
// Get normal to enclosing volume at boundary point in global coordinates

G4ThreeVector G4CMP::GetSurfaceNormal(const G4Step& step) {

  //RL: when we have nested geometries and the step is impinging upon an internal/daughter volume
  //from a parent volume, it seems that the returned normal is kind of wonky if one uses the pre-step volume
  //as a way to access that. So we need a bit of logic to handle these situations.

  //First, let's get our pre-point volume touchable and its associated solid  
  const G4VTouchable* preTouch = step.GetPreStepPoint()->GetTouchable();
  G4VSolid* preSolid = preTouch->GetVolume()->GetLogicalVolume()->GetSolid();
  G4ThreeVector pos = step.GetPostStepPoint()->GetPosition();

  //Rotate to local position -- use preTouch or postTouch? Leaving for now.
  G4cout << "REL in G4CMP::GetSurfaceNormal: preSolid name: " << preSolid->GetName() << G4endl;
  G4cout << "REL in G4CMP::GetSurfaceNormal: pre-rotation pos: " << pos.x() << ", " << pos.y() << ", " << pos.z() << G4endl;
  RotateToLocalPosition(preTouch, pos);
  G4cout << "REL in G4CMP::GetSurfaceNormal: post-rotation pos: " << pos.x() << ", " << pos.y() << ", " << pos.z() << G4endl;

  //Let's look at the distance to its exterior. This is done with the DistanceToOut function. If the distance
  //is nonzero, then we should be at a boundary with a daughter volume and we need to change the logic up
  //to use THAT boundary's surface normal
  G4cout << "REL in G4CMP::GetSurfaceNormal, distanceToOut: " << preSolid->DistanceToOut(pos) << G4endl;
  G4cout << "REL in G4CMP::GetSurfaceNormal, distanceToIn: " << preSolid->DistanceToIn(pos) << G4endl;
  double tolerance = 1.0e-9 * mm; //Need to not hardcode this -- find a better way to implement. 
  if( preSolid->DistanceToOut(pos) > tolerance ){
    G4cout << "REL in G4CMP::GetSurfaceNormal: On an internal (daughter) boundary of a volume." << G4endl;    
    const G4VTouchable* postTouch = step.GetPostStepPoint()->GetTouchable();   //REL added
    G4VSolid* postSolid = postTouch->GetVolume()->GetLogicalVolume()->GetSolid();   //REL added

    G4ThreeVector postSolidNorm = postSolid->SurfaceNormal(pos); //REL added
    G4cout << "REL in G4CMP::GetSurfaceNormal: postSolid norm, pre-rotation: " << postSolidNorm.x() << ", " << postSolidNorm.y() << ", " << postSolidNorm.z() << G4endl;
    RotateToGlobalDirection(postTouch,postSolidNorm); //REL added
    G4cout << "REL in G4CMP::GetSurfaceNormal: postSolid norm, post-rotation: " << postSolidNorm.x() << ", " << postSolidNorm.y() << ", " << postSolidNorm.z() << G4endl;
    return postSolidNorm;
  }
  else{
    G4cout << "REL in G4CMP::GetSurfaceNormal: On an external boundary of a volume." << G4endl;   
    G4ThreeVector norm = preSolid->SurfaceNormal(pos);
    G4cout << "REL in G4CMP::GetSurfaceNormal: norm, pre-rotation: " << norm.x() << ", " << norm.y() << ", " << norm.z() << G4endl;
    RotateToGlobalDirection(preTouch, norm);  
    G4cout << "REL in G4CMP::GetSurfaceNormal: norm, post-rotation: " << norm.x() << ", " << norm.y() << ", " << norm.z() << G4endl;
    return norm;
  }

}
*/

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


/*
//1. The touchable of the mother volume that the point is currently in
//2. The position of the point whose safety is to be computed
//3. The momentum direction of the last point (used for identifying what directions are "in-plane")
//4. A bool saying whether we're trying compute a safety from a boundary we're currently on. REL do we actually need this, or can we not calculate?
//5. A tangent vector indicating one angular bound for where to look
//6. A tangent vector indicating another angular bound for where to look
G4double G4CMP::ComputeConstrained2DSafety(const G4VTouchable* motherTouch, G4ThreeVector pos, G4ThreeVector momDir, bool safetyFromABoundary, G4ThreeVector tangVector1, G4ThreeVector tangVector2)
{
  //Pseudocode
  //0. Define output variables and make some copies
  G4double overallSafety = DBL_MAX;

  //Get the mother volume information
  G4VPhysicalVolume* motherPhys = motherTouch->GetVolume();
  G4LogicalVolume* motherLog = motherPhys->GetLogicalVolume();
  G4VSolid * motherSolid = motherLog->GetSolid();

  //Rotate the position and the momentum direction to the touchable's coordinates.
  RotateToLocalPosition(motherTouch, pos);
  RotateToLocalDirection(motherTouch, momDir);
  RotateToLocalDirection(motherTouch,tangVector1);
  RotateToLocalDirection(motherTouch,tangVector2);
  
  //1. First, get the shortest distance to the mother volume that we're in. ("DistanceToOut")
  G4double motherSafety = Compute2DSafetyInMotherVolume(motherTouch, motherPhys, motherLog, motherSolid, pos, safetyFromABoundary,tangVector1,tangVector2);
  if( motherSafety < overallSafety ){
    overallSafety = motherSafety;
  }

  //2. Next, loop through the daughter volumes in this mother volume and compute the distances to those.
  // Here, we don't actually use the tangVector1 and tangVector2 in the math, because of something that I think will be a reasonable
  //first approximation: if you're in a corner, the chances of a daughter boundary being on the same scale of how close you are to that
  //corner should be small. Some edge cases exist, where you can form a corner from a mother and a daugher boundary, in which case you
  //may land on the daughter, re-recognize that you're stuck, and then get re-ejected safely far from the corner. So at least in reasonably simple
  //geometries, this plus our mother safety should be fine. Need to rejigger this code to integrate it into the Get2D safety, since there are lots
  //of similarities/overlaps. (REL)
  for( int iD = 0; iD < motherLog->GetNoDaughters(); ++iD ){
    G4double daughterSafety = Compute2DSafetyToDaughterVolume(pos,momDir,motherLog,safetyFromABoundary,iD,tangVector1,tangVector2);
    if( daughterSafety < overallSafety ) overallSafety = daughterSafety;
  }
  
  //Return the safety
  return overallSafety; 
}

*/

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

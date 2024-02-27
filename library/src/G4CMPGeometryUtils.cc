/***********************************************************************\
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
// 20211001  Add utilities to get lattice from touchable, find valley close
//		to specified direction.

#include "G4CMPGeometryUtils.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPGlobalLocalTransformStore.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalVolume.hh"
#include "G4Navigator.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"
#include "G4VTouchable.hh"
#include "Randomize.hh"


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


// Get normal to enclosing volume at boundary point in global coordinates

G4ThreeVector G4CMP::GetSurfaceNormal(const G4Step& step) {
  const G4VTouchable* preTouch = step.GetPreStepPoint()->GetTouchable();
  G4VSolid* preSolid = preTouch->GetVolume()->GetLogicalVolume()->GetSolid();

  G4ThreeVector pos = step.GetPostStepPoint()->GetPosition();
  RotateToLocalPosition(preTouch, pos);

  G4ThreeVector norm = preSolid->SurfaceNormal(pos);
  RotateToGlobalDirection(preTouch, norm);
  
  return norm;
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

  G4VPhysicalVolume* pv = touch->GetVolume();
  G4LatticePhysical* lat = GetLattice(touch);
  if (!lat) {		// No lattice in touchable's volume, try pos instead
    pv = G4CMP::GetVolumeAtPoint(pos);
    lat = G4LatticeManager::GetLatticeManager()->GetLattice(pv);

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


// Get lattice associated with specific location in geometry

G4LatticePhysical* G4CMP::GetLattice(const G4VTouchable* touch) {
  if (!touch) return 0;

  G4VPhysicalVolume* pv = touch->GetVolume();
  return G4LatticeManager::GetLatticeManager()->GetLattice(pv);
}

// Find valley in crystal closest to specified direction
// NOTE:  Direction should be in GLOBAL coordinates, passed with touchable

G4int G4CMP::FindNearestValley(const G4VTouchable* touch, G4ThreeVector gdir) {
  if (!touch) return -1;

  RotateToLocalDirection(touch, gdir);
  return FindNearestValley(GetLattice(touch), gdir);
}

// Find valley in crystal closest to specified direction
// NOTE:  Direction should be in LOCAL coordinates, for use with lattice

G4int 
G4CMP::FindNearestValley(const G4LatticePhysical* lat, G4ThreeVector ldir) {
  if (!lat) return -1;

  ldir.setR(1.);		// Force to be unit vector
  lat->RotateToLattice(ldir);

  std::set<G4int> bestValley;	// Collect all best matches for later choice
  G4double align, bestAlign = -1.;
  G4ThreeVector pq = 0;
  for (size_t i=0; i<lat->NumberOfValleys(); i++) {
    pq = lat->MapPToP_Q(ldir, i);
    pq = pq.unit();
    align = fabs(lat->GetValleyAxis(i).dot(pq)); // Both unit vectors
    if (align > bestAlign) {
      bestValley.clear();
      bestAlign = align;
    }
    if (align >= bestAlign) bestValley.insert(i);
  }

  // Return best alignment, or pick from ambiguous choices
  G4int index = ( (bestValley.size() == 1) ? 0
		  : bestValley.size()*G4UniformRand() );

  return *std::next(bestValley.begin(), index);
}

/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef G4CMPGeometryUtils_hh
#define G4CMPGeometryUtils_hh 1

// $Id$
// File: G4CMPGeometryUtils.hh
//
// Description: Free standing helper functions for geometry based calculations.
//
// 20161107  Rob Agnese
// 20170605  Pass touchable from track, not just local PV
// 20170731  Add utility to get volume at (global) position
// 20170815  Add utility to shift position to avoid volume surfaces
// 20170913  Add utility to get electric field at (global) position
// 20170925  Add utility to create touchable at (global) position
// 20190226  Use local instance of G4Navigator to avoid corrupting tracking
// 20211001  Add utilities to get lattice from touchable, find valley close
//		to specified direction.
// 20240822  Add optional direction vector to GetVolumeAtPoint(), to break
//		ambiguities with surface points.

#include "G4ThreeVector.hh"

class G4LatticePhysical;
class G4Navigator;
class G4Step;
class G4Track;
class G4VPhysicalVolume;
class G4VTouchable;


namespace G4CMP {
  G4ThreeVector GetLocalDirection(const G4VTouchable* touch,
				  const G4ThreeVector& dir);
  
  G4ThreeVector GetLocalPosition(const G4VTouchable* touch,
				 const G4ThreeVector& pos);
  
  G4ThreeVector GetGlobalDirection(const G4VTouchable* touch,
				   const G4ThreeVector& dir);
  
  G4ThreeVector GetGlobalPosition(const G4VTouchable* touch,
				  const G4ThreeVector& pos);
  
  void RotateToLocalDirection(const G4VTouchable* touch,
			      G4ThreeVector& dir);
  
  void RotateToLocalPosition(const G4VTouchable* touch,
			     G4ThreeVector& pos);
  
  void RotateToGlobalDirection(const G4VTouchable* touch,
			       G4ThreeVector& dir);
  
  void RotateToGlobalPosition(const G4VTouchable* touch,
			      G4ThreeVector& pos);
  
  G4ThreeVector GetSurfaceNormal(const G4Step& step);

  G4Navigator* GetNavigator();		// Non-tracking for point finding

  // Find physical volume which contains specified point
  G4VPhysicalVolume* GetVolumeAtPoint(const G4ThreeVector& pos);

  // As above, but uses direction vector if point is on surface of volume(s)
  G4VPhysicalVolume* GetVolumeAtPoint(const G4ThreeVector& pos,
				      const G4ThreeVector& dir);

  // NOTE:  Transfers ownership to client
  G4VTouchable* CreateTouchableAtPoint(const G4ThreeVector& pos);

  G4ThreeVector ApplySurfaceClearance(const G4VTouchable* touch,
				      G4ThreeVector pos);

  G4LatticePhysical* GetLattice(const G4VTouchable* touch);

  // NOTE:  Direction should be in GLOBAL coordinates, passed with touchable
  G4int FindNearestValley(const G4VTouchable* touch, G4ThreeVector gdir);

  // NOTE:  Direction should be in LOCAL coordinates, for use with lattice
  G4int FindNearestValley(const G4LatticePhysical* lat, G4ThreeVector ldir);
}

#endif	/* G4CMPGeometryUtils_hh */

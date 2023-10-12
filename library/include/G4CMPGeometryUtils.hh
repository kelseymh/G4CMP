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

#include "G4ThreeVector.hh"

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

  G4VPhysicalVolume* GetVolumeAtPoint(const G4ThreeVector& pos);

  // NOTE:  Transfers ownership to client
  G4VTouchable* CreateTouchableAtPoint(const G4ThreeVector& pos);

  G4ThreeVector ApplySurfaceClearance(const G4VTouchable* touch,
				      G4ThreeVector pos);
}

#endif	/* G4CMPGeometryUtils_hh */

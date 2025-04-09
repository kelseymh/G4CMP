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
class G4StepPoint;
class G4Track;
class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VTouchable;
class G4VSolid;

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
  //G4ThreeVector GetSurfaceNormal(const G4StepPoint* stepPoint, const G4ThreeVector& guessedDirection );

  G4double Get2DSafety(const G4VTouchable* motherTouch,
		       G4ThreeVector pos,
		       G4ThreeVector momDir,
		       bool safetyFromABoundary,
		       bool forceSweepSafetyForDaughters = false,
		       G4ThreeVector surfaceNorm = G4ThreeVector(0,0,0),
		       G4ThreeVector tangVect1 = G4ThreeVector(0,0,0),
		       G4ThreeVector tangVect2 = G4ThreeVector(0,0,0));

  std::pair<G4double,G4ThreeVector> Get2DSafetyWithDirection(const G4VTouchable* motherTouch,
							     G4ThreeVector pos,
							     G4ThreeVector momDir,
							     bool safetyFromABoundary,
							     G4ThreeVector surfaceNorm = G4ThreeVector(0,0,0),
							     G4ThreeVector tangVect1 = G4ThreeVector(0,0,0),
							     G4ThreeVector tangVect2 = G4ThreeVector(0,0,0));

  
  G4double Compute2DSafetyToDaughterVolume(const G4ThreeVector & pos,
					   const G4ThreeVector & momDir,
					   G4LogicalVolume * motherLog,
					   bool safetyFromABoundary,
					   G4int daughterID,
					   G4ThreeVector & returnDir,
					   G4bool forceSweepSafetyForDaughters = false,
					   G4ThreeVector surfaceNorm = G4ThreeVector(0,0,0),
					   G4ThreeVector tangVect1 = G4ThreeVector(0,0,0),
					   G4ThreeVector tangVect2 = G4ThreeVector(0,0,0));
					   
  
  G4double Compute2DSafetyInMotherVolume(G4VSolid * motherSolid,
					 G4ThreeVector pos,
					 bool safetyFromABoundary,
					 G4ThreeVector & returnDir,
					 G4ThreeVector surfaceNorm = G4ThreeVector(0,0,0),
					 G4ThreeVector tangVect1 = G4ThreeVector(0,0,0),
					 G4ThreeVector tangVect2 = G4ThreeVector(0,0,0));
  
  G4double Compute2DSafetyFromABoundary(const G4VSolid * theVolSolid,
					G4ThreeVector pos,
					G4ThreeVector & returnDir,
					G4ThreeVector surfaceNorm,
					G4ThreeVector tangVect1,
					G4ThreeVector tangVect2,
					bool volIsMother);
  
  G4double GetSafetyInZ(const G4VTouchable* motherTouch,
			G4ThreeVector pos);
    
  
  G4double Compute2DMotherSafetyFromtheBulk(const G4VSolid * motherSolid,G4ThreeVector pos,G4ThreeVector & returnDir);

  G4double Compute2DDaughterSweptSafety(const G4VSolid* volDaughterSolid, G4ThreeVector pos,G4ThreeVector & returnDir);
  
  
  G4Navigator* GetNavigator();		// Non-tracking for point finding
  
  G4VPhysicalVolume* GetVolumeAtPoint(const G4ThreeVector& pos);

  // NOTE:  Transfers ownership to client
  G4VTouchable* CreateTouchableAtPoint(const G4ThreeVector& pos);

  G4ThreeVector ApplySurfaceClearance(const G4VTouchable* touch,
				      G4ThreeVector pos);

  G4double ComputeDotProductThreshold_Tang(int half_circle_nV);
  G4double ComputeDotProductThreshold_Norm(int half_circle_nV);
  
}

#endif	/* G4CMPGeometryUtils_hh */

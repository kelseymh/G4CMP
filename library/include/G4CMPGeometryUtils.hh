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
// 20250905  R. Linehan -- Add a function to robustify the random 2d vector
//           generation from 3D random vector
// 20251116  For G4 11, use #include "G4VTouchable.hh"

#include "G4ThreeVector.hh"
#include "G4VTouchable.hh"

class G4LatticePhysical;
class G4LogicalVolume;
class G4Navigator;
class G4Step;
class G4StepPoint;
class G4Track;
class G4VPhysicalVolume;
class G4VSolid;


namespace G4CMP {

  const G4ThreeVector nullVec(0,0,0); //For use in default initialization
  
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
					    
  G4ThreeVector GetSurfaceNormal(const G4Step& step,
				 const G4ThreeVector& incMomDir = nullVec );

  G4double Get2DSafety(const G4VTouchable* motherTouch,
                       const G4ThreeVector& pos_in,
                       const G4ThreeVector& momDir_in,
                       bool safetyFromABoundary,
                       bool forceSweepSafetyForDaughters = false,
                       const G4ThreeVector& surfaceNorm_in = nullVec,
                       const G4ThreeVector& tangVect1_in = nullVec,
                       const G4ThreeVector& tangVect2_in = nullVec);
 
  std::pair<G4double,G4ThreeVector>
  Get2DSafetyWithDirection(const G4VTouchable* motherTouch,
                           const G4ThreeVector& pos_in,
                           const G4ThreeVector& momDir_in,
                           bool safetyFromABoundary,
                           const G4ThreeVector& surfaceNorm_in = nullVec,
                           const G4ThreeVector& tangVect1_in = nullVec,
                           const G4ThreeVector& tangVect2_in = nullVec);  
  
  G4double
  Compute2DSafetyToDaughterVolume(const G4ThreeVector& pos,
                                  const G4ThreeVector& momDir,
                                  G4LogicalVolume * motherLog,
                                  bool safetyFromABoundary,
                                  G4int daughterID,
                                  G4ThreeVector & returnDir,
                                  G4bool forceSweepSafetyForDaughters = false,
                                  const G4ThreeVector& surfaceNorm_in = nullVec,
                                  const G4ThreeVector& tangVect1_in = nullVec,
                                  const G4ThreeVector& tangVect2_in = nullVec );

  G4double
  Compute2DSafetyInMotherVolume(G4VSolid * motherSolid,
                                const G4ThreeVector & pos,
                                bool safetyFromABoundary,
                                G4ThreeVector & returnDir,
                                const G4ThreeVector& surfaceNorm_in = nullVec,
                                const G4ThreeVector& tangVect1_in = nullVec,
                                const G4ThreeVector& tangVect2_in = nullVec );

  G4double 
  Compute2DSafetyFromABoundary(const G4VSolid * theVolSolid,
                               const G4ThreeVector & pos,
                               G4ThreeVector & returnDir,
                               const G4ThreeVector& surfaceNorm_in,
                               const G4ThreeVector& tangVect1_in,
                               const G4ThreeVector& tangVect2_in,
                               bool volIsMother);
  
  G4double GetSafetyInZ(const G4VTouchable* motherTouch,
                        const G4ThreeVector& pos_in);
  
  G4double Compute2DMotherSafetyFromtheBulk(const G4VSolid * motherSolid,
                                            const G4ThreeVector & pos,
                                            G4ThreeVector & returnDir);

  G4double Compute2DDaughterSweptSafety(const G4VSolid* volDaughterSolid,
                                        const G4ThreeVector& pos_in,
                                        G4ThreeVector & returnDir);
  
  G4Navigator* GetNavigator();		// Non-tracking for point finding
  
  G4VPhysicalVolume* GetVolumeAtPoint(const G4ThreeVector& pos);

  // NOTE:  Transfers ownership to client
  G4VTouchable* CreateTouchableAtPoint(const G4ThreeVector& pos);

  G4ThreeVector ApplySurfaceClearance(const G4VTouchable* touch,
				      G4ThreeVector pos);


  G4double ComputeDotProductThreshold_Tang(int half_circle_nV);
  G4double ComputeDotProductThreshold_Norm(int half_circle_nV);
  G4LatticePhysical* GetLattice(const G4VTouchable* touch);

  // NOTE:  Direction should be in GLOBAL coordinates, passed with touchable
  G4int FindNearestValley(const G4VTouchable* touch, G4ThreeVector gdir);

  // NOTE:  Direction should be in LOCAL coordinates, for use with lattice
  G4int FindNearestValley(const G4LatticePhysical* lat, G4ThreeVector ldir);

  G4ThreeVector RobustifyRandomDirIn2D(G4ThreeVector returnDir);  
}

#endif	/* G4CMPGeometryUtils_hh */

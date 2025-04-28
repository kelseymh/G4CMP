/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPSolidUtils.hh
/// \brief This class contains custom utilities for G4VSolids that are not
/// directly included in the G4VSolid class. Some example of utilites are:
/// 1) Adjusting a position to the nearest surface point without the surface normal.
/// 2) Adjusting a particle to the nearest edge along the surface.
///
/// The user will need to create this class directly with a G4VSolid and/or
/// touchable, or can initialize these objects later on.
///
/// Note: All coordinates/vectors passed into functions must be in the solid's local coordinate system.
//
// $Id$
//
// 20250424  G4CMP-465 -- Create G4CMPSolidUtils class.

#ifndef G4CMPSolidUtils_hh
#define G4CMPSolidUtils_hh 1

#include "globals.hh"
#include "G4AffineTransform.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

class G4VSolid;
class G4VTouchable;

class G4CMPSolidUtils {
  public:
    // Default constructor
    G4CMPSolidUtils() : theSolid(0), theTransform(G4AffineTransform()), verboseLevel(0),
                        verboseLabel("G4CMPSolidUtils") {;}

    // Direct constructor with solid & transform for client code in local frame
    G4CMPSolidUtils(const G4VSolid* solid, G4int verbose, G4String vLabel);
    G4CMPSolidUtils(const G4VSolid* solid, const G4AffineTransform& trans,
                    G4int verbose, G4String vLabel);
    G4CMPSolidUtils(const G4VSolid* solid, const G4RotationMatrix& rot,
                    const G4ThreeVector& disp, G4int verbose, G4String vLabel);

    // Direct constructor with touchable for client code in global frame
    G4CMPSolidUtils(const G4VTouchable* touch, G4int verbose, G4String vLabel);

    // Copy operation
    G4CMPSolidUtils& operator=(const G4CMPSolidUtils& right);

    // Set and get member variables
    void SetSolid(const G4VSolid* solid) {
      theSolid = solid;
    }

    void SetTransform(const G4AffineTransform& trans) {
      theTransform = trans;
    }

    void SetTransform(const G4RotationMatrix& rot,
                      const G4ThreeVector& disp) {
      theTransform.SetNetRotation(rot);
      theTransform.SetNetTranslation(disp);
    }

    void SetVerboseLevel(G4int verbose) {
      verboseLevel = verbose;
    }

    void SetVerboseLevel(G4int verbose, G4String vLabel) {
      verboseLevel = verbose;
      verboseLabel = vLabel;
    }

    void SetVerboseLabel(G4String vLabel) {
      verboseLabel = vLabel;
    }

    const G4VSolid* GetSolid() const { return theSolid; }
    const G4AffineTransform GetTransform() const { return theTransform; }
    G4int GetVerboseLevel() const { return verboseLevel; }
    G4String GetVerboseLabel() const {return verboseLabel; }

    // Get the distance to solid object surface with and without direction
    G4double GetDistanceToSolid(const G4ThreeVector& localPos) const;
    G4double GetDistanceToSolid(const G4ThreeVector& localPos,
                                const G4ThreeVector& dir) const;

    // Get the direction for minimum distance to the solid object's surface
    // localPos and dir need to be in the solid's local coordinate system
    G4ThreeVector GetDirectionToSolid(const G4ThreeVector& localPos) const;

    // Modify dir in place to direction of minimum distance to the solid object
    // localPos and dir need to be in the solid's local coordinate system
    void RotateDirectionToSolid(const G4ThreeVector& localPos,
                                G4ThreeVector& dir) const;

    // Efficiently find direction with min distance to surface
    // localPos must be in the solid's local coordinate system
    void OptimizeSurfaceAdjustAngle(const G4ThreeVector& localPos,
                                    G4double& theta0, G4double& phi0,
                                    const G4int angOption,
                                    const G4double minDist) const;

    // Modifies localPos in place to closest surface position
    // localPos and dir must be in the solid's local coordinate system
    void AdjustToClosestSurfacePoint(G4ThreeVector& localPos) const;
    void AdjustToClosestSurfacePoint(G4ThreeVector& localPos,
                                     const G4ThreeVector& dir) const;

    // localPos and dir must be in the solid's local coordinate system
    G4ThreeVector GetClosestSurfacePoint(const G4ThreeVector& localPos) const;
    G4ThreeVector GetClosestSurfacePoint(const G4ThreeVector& localPos,
                                         const G4ThreeVector& dir) const;

    // Modifies localPos in place to the nearest edge along vTan
    // localPos and vTan must be in the solid's local coordinate system
    void AdjustToEdgePosition(const G4ThreeVector& vTan,
                              G4ThreeVector& localPos,
                              G4double high, const G4int curvedSurf) const;

    // localPos and vTan must be in the solid's local coordinate system
    G4ThreeVector GetEdgePosition(const G4ThreeVector& vTan,
                                  const G4ThreeVector& localPos,
                                  G4double high, const G4int curvedSurf) const;

    // Modifies vTan and surfNorm in place
    // localPos, vTan, and surfNorm must be in the solid's local coordinate system
    void ReflectAgainstEdge(G4ThreeVector& vTan,
                            const G4ThreeVector& localPos,
                            G4ThreeVector& surfNorm) const;

  protected:
    void TransformLocalToGlobal(G4ThreeVector& point, G4ThreeVector& dir) const;
    void TransformGlobalToLocal(G4ThreeVector& point, G4ThreeVector& dir) const;
    void TransformToGlobalPoint(G4ThreeVector& point) const;
    void TransformToLocalPoint(G4ThreeVector& point) const;
    void TransformToGlobalDirection(G4ThreeVector& dir) const;
    void TransformToLocalDirection(G4ThreeVector& dir) const;

    const G4ThreeVector GetLocalPosition(const G4ThreeVector& pos) const {
      G4ThreeVector localPos = pos;
      TransformToLocalPoint(localPos);
      return localPos;
    }

    const G4ThreeVector GetGlobalPosition(const G4ThreeVector& localPos) {
      G4ThreeVector globalPos = localPos;
      TransformToGlobalPoint(globalPos);
      return globalPos;
    }

    const G4ThreeVector GetLocalDirection(const G4ThreeVector& dir) const {
      G4ThreeVector localDir = dir;
      TransformToLocalDirection(localDir);
      return localDir;
    }

    const G4ThreeVector GetGlobalDirection(const G4ThreeVector& localDir) const {
      G4ThreeVector globalDir = localDir;
      TransformToGlobalDirection(globalDir);
      return globalDir;
    }

  private:
    const G4VSolid* theSolid;
    G4AffineTransform theTransform;
    G4int verboseLevel;
    G4String verboseLabel;
};

#endif	/* G4CMPSolidUtils_hh */
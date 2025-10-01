/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4PhononReflection.cc
/// \brief This class contains custom utilities for G4VSolids that are not
/// directly included in the G4VSolid class. Some example of utilites are:
/// 1) Adjusting a position to the nearest surface point without the surface
///    normal.
/// 2) Adjusting a particle to the nearest edge along the surface.
/// 3) Reflecting a tangent vector against a "hard" edge.
///
/// The user will need to create this class directly with a G4VSolid (and
/// transform or a touchable). If the transform or touchable is not
/// provided, an indentity transform is used.
///
/// Note: All coordinates/vectors passed into functions must be in the global
/// coordinate system. If a global-to-local transform is supplied, the input
/// coordinates and vectors, and any output results, will be in the coordinate
/// system of that transform.
///
//
// $Id$
//
// 20250424  G4CMP-465 -- Create G4CMPSolidUtils class.
// 20250429  G4CMP-461 -- Add function for skipping detector flats.
// 20250430  N. Tenpas -- Add function for getting distance to bounding box.

#include "G4CMPSolidUtils.hh"
#include "G4AffineTransform.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4VSolid.hh"
#include "G4VTouchable.hh"
#include "G4UnitsTable.hh"


// Direct constructors

G4CMPSolidUtils::G4CMPSolidUtils(const G4VSolid* solid,
                                 G4int verbose, const G4String& vLabel)
  : theSolid(solid), theTransform(G4AffineTransform()), verboseLevel(verbose),
    verboseLabel(vLabel) {;}

G4CMPSolidUtils::G4CMPSolidUtils(const G4VSolid* solid,
                                 const G4AffineTransform& trans,
                                 G4int verbose, const G4String& vLabel)
  : theSolid(solid), theTransform(trans), verboseLevel(verbose),
    verboseLabel(vLabel) {;}

G4CMPSolidUtils::G4CMPSolidUtils(const G4VSolid* solid,
                                 const G4RotationMatrix& rot,
                                 const G4ThreeVector& disp,
                                 G4int verbose, const G4String& vLabel)
  : theSolid(solid), theTransform(G4AffineTransform(rot, disp)),
    verboseLevel(verbose), verboseLabel(vLabel) {;}

G4CMPSolidUtils::G4CMPSolidUtils(const G4VTouchable* touch, G4int verbose,
                                 const G4String& vLabel)
  : theSolid(touch->GetSolid()),
    theTransform(G4AffineTransform(touch->GetRotation(),
				   touch->GetTranslation())),
    verboseLevel(verbose), verboseLabel(vLabel) {;}


// Copy operation

G4CMPSolidUtils&
G4CMPSolidUtils::operator=(const G4CMPSolidUtils& right) {
  theSolid = right.theSolid;
  theTransform = right.theTransform;
  verboseLevel = right.verboseLevel;
  verboseLabel = right.verboseLabel;
  return *this;
}


// Get the shortest distance to solid

G4double G4CMPSolidUtils::GetDistanceToSolid(const G4ThreeVector& pos) const {
  G4ThreeVector localPos = GetLocalPosition(pos);
  return theSolid->DistanceToOut(localPos) + theSolid->DistanceToIn(localPos);
}

// Get distance to solid along direction

G4double G4CMPSolidUtils::GetDistanceToSolid(const G4ThreeVector& pos,
                                             const G4ThreeVector& dir) const {
  G4ThreeVector localPos = GetLocalPosition(pos);
  G4ThreeVector localDir = GetLocalDirection(dir);

  if (theSolid->Inside(localPos) == kInside) {
    return theSolid->DistanceToOut(localPos, localDir);
  } else {
    return theSolid->DistanceToIn(localPos, localDir);
  }
}


// Get the direction for the shortest distance to a solid

G4ThreeVector G4CMPSolidUtils::GetDirectionToSolid(const G4ThreeVector& pos) const {
  G4ThreeVector optDir(0,0,0);
  RotateDirectionToSolid(pos, optDir);
  return optDir;
}

// Adjust dir in place to direction with shortest distance to a solid

void G4CMPSolidUtils::RotateDirectionToSolid(const G4ThreeVector& pos,
                                             G4ThreeVector& dir) const {
  // Angles to be adjusted in place by OptimizeSurfaceAdjustAngle
  // Start at theta = pi/2 so "fit" can be determined by phi
  G4double bestTheta = pi / 2;
  G4double bestPhi = 0;

  G4double minDist = GetDistanceToSolid(pos);

  // Try to optimize phi first
  OptimizeSurfaceAdjustAngle(pos, bestTheta, bestPhi, 1, minDist);
  OptimizeSurfaceAdjustAngle(pos, bestTheta, bestPhi, 0, minDist);

  dir.setRThetaPhi(minDist, bestTheta, bestPhi);

  // Return a null vector if no direction was found
  if (theSolid->Inside(GetLocalPosition(pos + dir)) != kSurface) dir.set(0,0,0);
  else dir.setMag(1.);
}


// Find the direction for minimum distance to surface
// Adjust theta0 and phi0 in place with best match
// Everything remains in global coordinates
void G4CMPSolidUtils::OptimizeSurfaceAdjustAngle(const G4ThreeVector& pos,
                                                 G4double& theta0,
                                                 G4double& phi0,
                                                 const G4int angOption,
                                                 const G4double minDist) const {
  // Constants used for Golden Section searches
  G4double const tolerance = 1e-12*rad;
  G4double const gRatio = (std::sqrt(5) + 1) / 2;

  // Define search quantities for finding optimal angle
  G4double a = 0;
  G4double b = (angOption == 0) ? pi : twopi;
  G4double x1 = a + (b - a) / gRatio;
  G4double x2 = b - (b - a) / gRatio;
  G4ThreeVector dir1(0,0,0);
  G4ThreeVector dir2(0,0,0);

  // Initial search directions
  if (angOption == 0) {
    // Optimize theta
    dir1.setRThetaPhi(minDist, x1, phi0);
    dir2.setRThetaPhi(minDist, x2, phi0);
  } else {
    // Optimize phi
    dir1.setRThetaPhi(minDist, theta0, x1);
    dir2.setRThetaPhi(minDist, theta0, x2);
  }

  // Distances for each direction
  G4double dist1 = GetDistanceToSolid(pos + dir1);
  G4double dist2 = GetDistanceToSolid(pos + dir2);

  // Golden section search for optimizing angle
  while (b - a > tolerance) {
    // Adjust bounds
    if (dist1 < dist2) {
      // Shift up
      a = x2;
      x2 = x1;
      dist2 = dist1;
      x1 = a + (b - a) / gRatio;
      if (angOption == 0) {
        // Optimize theta
        dir1.setRThetaPhi(minDist, x1, phi0);
      } else {
        // Optimize phi
        dir1.setRThetaPhi(minDist, theta0, x1);
      }
      dist1 = GetDistanceToSolid(pos + dir1);
    } else {
      // Shift down
      b = x1;
      x1 = x2;
      dist1 = dist2;
      x2 = b - (b - a) / gRatio;
      if (angOption == 0) {
        // Optimize theta
        dir2.setRThetaPhi(minDist, x2, phi0);
      } else {
        // Optimize phi
        dir2.setRThetaPhi(minDist, theta0, x2);
      }
      dist2 = GetDistanceToSolid(pos + dir2);
    }
  }
  G4double optimalAng = (x1 + x2) / 2;

  // Adjust angle in place
  if (angOption == 0) { theta0 = optimalAng; }
  else { phi0 = optimalAng; }
}


// Get closest position on surface of solid

G4ThreeVector G4CMPSolidUtils::GetClosestSurfacePoint(const G4ThreeVector& pos) const {
  G4ThreeVector surfacePoint = pos;
  AdjustToClosestSurfacePoint(surfacePoint);
  return surfacePoint;
}

// Adjust to surface position closest to pos without direction

void G4CMPSolidUtils::
AdjustToClosestSurfacePoint(G4ThreeVector& pos) const {
  // Do nothing if already on surface
  if (theSolid->Inside(GetLocalPosition(pos)) == kSurface) return;

  G4double minDist = GetDistanceToSolid(pos);
  G4ThreeVector optDir = minDist * GetDirectionToSolid(pos);

  // Only return valid positions on surface
  if (theSolid->Inside(GetLocalPosition(pos + optDir)) == kSurface) {
    pos += optDir;
  } else {
    if (verboseLevel>2) {
      G4cout << verboseLabel << "::AdjustToClosestSurfacePoint"
	     << ": Surface point not found from initial position "
	     << pos << G4endl;
    }
    pos.set(kInfinity,kInfinity,kInfinity);
  }
}

// Get closest position on surface of solid along direction

G4ThreeVector G4CMPSolidUtils::GetClosestSurfacePoint(const G4ThreeVector& pos,
                                                      const G4ThreeVector& dir) const {
G4ThreeVector surfacePoint = pos;
AdjustToClosestSurfacePoint(surfacePoint, dir);
return surfacePoint;
}

// Adjust to surface position closest to pos using a provided direction

void G4CMPSolidUtils::AdjustToClosestSurfacePoint(G4ThreeVector& pos,
                                                  const G4ThreeVector& dir) const {
  // Do nothing if on surface
  if (theSolid->Inside(GetLocalPosition(pos)) == kSurface) return;

  // Adjust position in place along provided direction
  G4double surfAdjust = GetDistanceToSolid(pos, dir);
  pos += surfAdjust * dir;

  if (theSolid->Inside(GetLocalPosition(pos)) != kSurface && verboseLevel>2) {
    G4cout << verboseLabel << "::AdjustToClosestSurfacePoint:"
	   << " WARNING: No surface found along " << dir
	   << ". Adjustment length: " << G4BestUnit(surfAdjust, "Length")
	   << G4endl;
  }
}


// Get edge position along vTan from pos
G4ThreeVector G4CMPSolidUtils::GetEdgePosition(const G4ThreeVector& vTan,
                                               const G4ThreeVector& pos,
                                               G4double high,
                                               const G4int curvedSurf) const {
  G4ThreeVector edgePosition = pos;
  AdjustToEdgePosition(vTan, edgePosition, high, curvedSurf);
  return edgePosition;
}

// Do a binary search to find the closest point toward the edge along vTan
// Modifies pos in place

void G4CMPSolidUtils::
AdjustToEdgePosition(const G4ThreeVector& vTan,
                     G4ThreeVector& pos, G4double high,
                     const G4int curvedSurf) const {
  EInside isIn = theSolid->Inside(GetLocalPosition(pos));
  G4double low = 0.0*um;
  G4double mid = 0;
  G4ThreeVector originalPos = pos;
  G4double tolerance = theSolid->GetTolerance();

  // Binary search to bring surface point to edge
  G4int maxItr = 100;
  G4int itr = 0;
  while ((high - low > tolerance || isIn != kSurface) && ++itr < maxItr) {
    mid = 0.5 * (low + high);
    pos = originalPos + mid * vTan;

    // Modify pos in place to surface if adjusting on curved surfaces
    if (curvedSurf == 1) AdjustToClosestSurfacePoint(pos);
    isIn = theSolid->Inside(GetLocalPosition(pos));

    if (isIn == kSurface) low = mid; // Move out
    else high = mid; // Move in
  }

  if (verboseLevel>2) {
    G4cout << verboseLabel << "::AdjustToEdgePosition"
	   << ": initialPos = " << originalPos
	   << ", vTan = " << vTan
	   << ", finalPos = " << pos << G4endl;
  }
}


// Reflect vTan against an edge
// Modifies vTan and surfNorm in place

void G4CMPSolidUtils::
ReflectAgainstEdge(G4ThreeVector& vTan, const G4ThreeVector& pos,
                   G4ThreeVector& surfNorm) const {
  G4ThreeVector localPos = GetLocalPosition(pos);
  TransformToLocalDirection(vTan);

  // Get normal of both surfaces at this edge point by stepping with normAdjust
  G4double normAdjust = 1*nm;
  G4double vTanMag = vTan.mag();
  G4ThreeVector norm1 = theSolid->SurfaceNormal(pos);
  G4ThreeVector norm2 = theSolid->SurfaceNormal(pos - normAdjust*norm1);
  // Try to fix norm1 if we didn't get the normals for a proper reflection
  norm1 = (norm1 * norm2 > 0) ? theSolid->SurfaceNormal(pos - normAdjust*norm2) : norm1;

  G4ThreeVector edgeVec(0,0,0);
  G4ThreeVector refNorm(0,0,0);

  // Only do reflection if at an edge
  if (norm1 * norm2 <= 0) {
    surfNorm = (surfNorm*norm1 > surfNorm*norm2) ? norm1 : norm2;
    // Project vTan to be orthogonal to one normal
    (vTan -= surfNorm * (vTan * surfNorm)).setMag(vTanMag);

    // Get the edge vector
    edgeVec = norm1.cross(norm2).unit();

    // Find the normal for the reflection "surface"
    refNorm = (edgeVec.cross(surfNorm)).unit();
    if (refNorm * vTan < 0) refNorm *= -1;

    // Reflect vTan against reflection "surface"
    vTan -= 2*((vTan * refNorm) * refNorm);
  }

  if (verboseLevel>2) {
    G4cout << verboseLabel << "::ReflectAgainstEdge"
	   << ": pos = " << pos
	   << ", vTan_0 = " << vTan - 2*((vTan * -refNorm) * -refNorm)
	   << ", edgeVector = " << edgeVec
	   << ", refNorm = " << refNorm
	   << ", norm1 = " << norm1
	   << ", norm2 = " << norm2
	   << ", vTan_f = " << vTan << G4endl;
  }

  TransformToGlobalDirection(vTan);
}


// Skip flats on surface (where the normal remains the same)
void G4CMPSolidUtils::AdjustOffFlats(G4ThreeVector& pos, G4ThreeVector& vTan,
                                     const G4double flatStepSize,
                                     G4ThreeVector& surfNorm, G4int count) {
  // Debugging variables
  G4ThreeVector originalPos = pos;
  G4ThreeVector originalV = vTan;

  // We must be under the consecutive call limit
  const G4int limit = 50;
  if (count > limit) {
    pos.set(kInfinity, kInfinity, kInfinity);
    return;
  }

  // We must start in the surface
  if (theSolid->Inside(GetLocalPosition(pos)) != kSurface) return;

  // Project vTan to guarantee it's tangent to the surface
  if (vTan * surfNorm != 0) {
    G4double vTanMag = vTan.mag();
    (vTan -= surfNorm * (vTan * surfNorm)).setMag(vTanMag);

    if (verboseLevel>2) {
      G4cout << verboseLabel << "::AdjustOffFlats"
	     << ": Tangent vector: " << originalV
	     << " is not orthogonal to surface normal: " << surfNorm
	     << ". New Tangent vector: " << vTan << G4endl;
    }
  }

  // Adjusts pos in place to the edge of the flat
  AdjustToEdgePosition(vTan.unit(), pos, flatStepSize, 0);

  // Test whether the next step is at an edge or on the curved wall
  const G4double surfAdjust = GetDistanceToSolid(pos + 1*nm*vTan.unit(), -surfNorm);
  G4ThreeVector localTrial = GetLocalPosition(pos + 1*nm*vTan.unit());
  // Adjust to surface
  localTrial -= surfAdjust * GetLocalDirection(surfNorm);

  if (verboseLevel>2) {
    G4cout << verboseLabel << "::AdjustOffFlats"
	   << ": Original Position = " << originalPos
	   << ", Final Position = " << pos
	   << ", Step Direction = " << vTan
	   << ", Surface Normal = " << surfNorm
	   << ", Surface Adjustment at pos = " << surfAdjust
	   << ", Skipper Step Size = " << flatStepSize / mm << " mm"
	   << ", Consecutive iterations count: " << count << G4endl;
  }

  // At a hard edge - reflect kTan and repeat
  if (theSolid->Inside(localTrial) != kSurface) {
    ReflectAgainstEdge(vTan, pos, surfNorm);
    AdjustOffFlats(pos, vTan, flatStepSize, surfNorm, ++count);
  }
}


// Find the boundary box position

G4double G4CMPSolidUtils::GetDistToBB(const G4ThreeVector pos,
                                      const G4ThreeVector vTan) const {
  G4ThreeVector bbMin, bbMax;
  theSolid->BoundingLimits(bbMin, bbMax);

  G4double dist_to_bb = 0.0;
  G4ThreeVector bbPos = pos;

  if (fabs(vTan.z()) > 1e-7) {
    // calculate distance to zMin
    dist_to_bb = (bbMin.z() - pos.z()) / vTan.z();

    if (dist_to_bb < 0) {
      // calculate distance to zMax
      dist_to_bb = (bbMax.z() - pos.z()) / vTan.z();
    }
  }

  if (fabs(vTan.x()) > 1e-7) {
    // do checks for Y-flats
    if (bbPos.x() < bbMin.x()) {
      // calculate distance to xMin
      dist_to_bb = (pos.x() - bbMin.x()) / vTan.x();
    } else if (bbPos.x() > bbMax.x()) {
      // calculate distance to xMax
      dist_to_bb = (bbMax.x() - pos.x()) / vTan.x();
    }
  }

  if (fabs(vTan.y()) > 1e-7) {
    // do checks for X-flats
    if (bbPos.y() < bbMin.y()) {
      // calculate distance to yMin
      dist_to_bb = (pos.y() - bbMin.y()) / vTan.y();
    } else if (bbPos.y() > bbMax.y()) {
      // calculate distance to yMax
      dist_to_bb = (bbMax.y() - pos.y()) / vTan.y();
    }
  }

  return dist_to_bb;
}


// Coordinate transformations

// Both position and direction transforms
void G4CMPSolidUtils::TransformGlobalToLocal(G4ThreeVector& point,
                                             G4ThreeVector& dir) const {
  TransformToLocalPoint(point);
  TransformToLocalDirection(dir);
}

void G4CMPSolidUtils::TransformLocalToGlobal(G4ThreeVector& point,
                                             G4ThreeVector& dir) const {
  TransformToGlobalPoint(point);
  TransformToGlobalDirection(dir);
}

// Position transforms
void G4CMPSolidUtils::TransformToGlobalPoint(G4ThreeVector& point) const {
  theTransform.ApplyPointTransform(point);
}

void G4CMPSolidUtils::TransformToLocalPoint(G4ThreeVector& point) const {
  theTransform.Inverse().ApplyPointTransform(point);
}

// Direction transforms
void G4CMPSolidUtils::TransformToGlobalDirection(G4ThreeVector& dir) const {
  theTransform.ApplyAxisTransform(dir);
}

void G4CMPSolidUtils::TransformToLocalDirection(G4ThreeVector& dir) const{
  theTransform.Inverse().ApplyAxisTransform(dir);
}

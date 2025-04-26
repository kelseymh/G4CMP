/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4PhononReflection.cc
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

#include "G4CMPSolidUtils.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4VSolid.hh"
#include "G4UnitsTable.hh"


// Copy operations

G4CMPSolidUtils::G4CMPSolidUtils(const G4VSolid* solid)
  : theSolid(solid), verboseLevel(0) {;}

G4CMPSolidUtils&
G4CMPSolidUtils::operator=(const G4CMPSolidUtils& right) {
  theSolid = right.theSolid;
  verboseLevel = right.verboseLevel;
  return *this;
}

// Initialize the solid / solid & verboseLevel

void G4CMPSolidUtils::Initialize(const G4VSolid* solid) {
  theSolid = solid;
}

void G4CMPSolidUtils::Initialize(const G4VSolid* solid, G4int verbose) {
  theSolid = solid;
  SetVerboseLevel(verbose);
}


// Get the shortest distance to solid

G4double G4CMPSolidUtils::GetDistanceToSolid(const G4ThreeVector& localPos) const {
  if (theSolid->Inside(localPos) == kInside) return theSolid->DistanceToOut(localPos);
  return theSolid->DistanceToIn(localPos);
}

// Get distance to solid along direction

G4double G4CMPSolidUtils::GetDistanceToSolid(const G4ThreeVector& localPos,
                                             const G4ThreeVector& dir) const {
  if (theSolid->Inside(localPos) == kInside) return theSolid->DistanceToOut(localPos, dir);
  return theSolid->DistanceToIn(localPos, dir);
}


// Get the direction for the shortest distance to a solid

G4ThreeVector G4CMPSolidUtils::GetDirectionToSolid(const G4ThreeVector& localPos) const {
  G4ThreeVector optDir(0,0,0);
  RotateDirectionToSolid(localPos, optDir);
  return optDir;
}

// Adjust dir in place to direction with shortest distance to a solid

void G4CMPSolidUtils::RotateDirectionToSolid(const G4ThreeVector& localPos,
                                             G4ThreeVector& dir) const {
  // Angles to be adjusted in place by OptimizeSurfaceAdjustAngle
  // Start at theta = pi/2 so "fit" can be determined by phi
  G4double bestTheta = pi / 2;
  G4double bestPhi = 0;

  G4double minDist = GetDistanceToSolid(localPos);

  // Try to optimize phi first
  OptimizeSurfaceAdjustAngle(localPos, bestTheta, bestPhi, 1, minDist);
  OptimizeSurfaceAdjustAngle(localPos, bestTheta, bestPhi, 0, minDist);

  dir.set(minDist*sin(bestTheta)*cos(bestPhi),
          minDist*sin(bestTheta)*sin(bestPhi),
          minDist*cos(bestTheta));

  // Return a null vector if no direction was found
  if (theSolid->Inside(localPos + dir) != kSurface) dir.set(0,0,0);
  else dir.setMag(1.);
}


// Find the direction for minimum distance to surface
// Adjust theta0 and phi0 in place with best match

void G4CMPSolidUtils::OptimizeSurfaceAdjustAngle(const G4ThreeVector& localPos,
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
  G4double coeffX = 0, coeffY = 0, coeffZ = 0;

  // Initial search directions
  if (angOption == 0) {
    // Optimize theta
    coeffX = minDist*cos(phi0);
    coeffY = minDist*sin(phi0);
    coeffZ = minDist;
    dir1.set(coeffX*sin(x1), coeffY*sin(x1), coeffZ*cos(x1));
    dir2.set(coeffX*sin(x2), coeffY*sin(x2), coeffZ*cos(x2));
  } else {
    // Optimize phi
    coeffX = minDist*sin(theta0);
    coeffY = minDist*sin(theta0);
    coeffZ = minDist*cos(theta0);
    dir1.set(coeffX*cos(x1), coeffY*sin(x1), coeffZ);
    dir2.set(coeffX*cos(x2), coeffY*sin(x2), coeffZ);
  }

  // Distances for each direction
  G4double dist1 = GetDistanceToSolid(localPos + dir1);
  G4double dist2 = GetDistanceToSolid(localPos + dir2);

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
        dir1.set(coeffX*sin(x1), coeffY*sin(x1), coeffZ*cos(x1));
      } else {
        // Optimize phi
        dir1.set(coeffX*cos(x1), coeffY*sin(x1), coeffZ);
      }
      dist1 = GetDistanceToSolid(localPos + dir1);
    } else {
      // Shift down
      b = x1;
      x1 = x2;
      dist1 = dist2;
      x2 = b - (b - a) / gRatio;
      if (angOption == 0) {
        // Optimize theta
        dir2.set(coeffX*sin(x2), coeffY*sin(x2), coeffZ*cos(x2));
      } else {
        // Optimize phi
        dir2.set(coeffX*cos(x2), coeffY*sin(x2), coeffZ);
      }
      dist2 = GetDistanceToSolid(localPos + dir2);
    }
  }
  G4double optimalAng = (x1 + x2) / 2;

  // Adjust angle in place
  if (angOption == 0) { theta0 = optimalAng; }
  else { phi0 = optimalAng; }
}


// Adjust to surface position closest to localPos without direction

void G4CMPSolidUtils::
AdjustToClosestSurfacePoint(G4ThreeVector& localPos) const {
  // Do nothing if already on surface
  if (theSolid->Inside(localPos) == kSurface) return;

  G4double minDist = GetDistanceToSolid(localPos);
  G4ThreeVector optDir = minDist * GetDirectionToSolid(localPos);

  // Only return valid positions on surface
  if (theSolid->Inside(localPos + optDir) == kSurface) {
    localPos += optDir;
  } else {
    localPos.set(kInfinity,kInfinity,kInfinity);
  }
}

// Get closest position on surface of solid

G4ThreeVector G4CMPSolidUtils::GetClosestSurfacePoint(const G4ThreeVector& localPos) const {
  G4ThreeVector surfacePoint = localPos;
  AdjustToClosestSurfacePoint(surfacePoint);
  return surfacePoint;
}

// Adjust to surface position closest to localPos using a provided direction

void G4CMPSolidUtils::AdjustToClosestSurfacePoint(G4ThreeVector& localPos,
                                                  const G4ThreeVector& dir) const {
  // Do nothing if on surface
  if (theSolid->Inside(localPos) == kSurface) return;

  // Adjust position in place along provided direction
  G4double surfAdjust = GetDistanceToSolid(localPos, dir);
  localPos += surfAdjust * dir;

  if (theSolid->Inside(localPos) != kSurface && verboseLevel) {
    G4cout << "WARNING: No surface found along " << dir
    << ". Adjustment length: " << G4BestUnit(surfAdjust, "Length") << G4endl;
  }
}

// Get closest position on surface of solid along direction

G4ThreeVector G4CMPSolidUtils::GetClosestSurfacePoint(const G4ThreeVector& localPos,
                                                      const G4ThreeVector& dir) const {
  G4ThreeVector surfacePoint = localPos;
  AdjustToClosestSurfacePoint(surfacePoint, dir);
  return surfacePoint;
}


// Do a binary search to find the closest point toward the edge along vTan
// Modifies localPos in place

void G4CMPSolidUtils::
AdjustToEdgePosition(const G4ThreeVector& vTan,
                     G4ThreeVector& localPos, G4double high,
                     const G4int curvedSurf) const {
  EInside isIn = theSolid->Inside(localPos);
  G4double low = 0.0*um;
  G4double mid = 0;
  G4ThreeVector originalPos = localPos;
  G4double tolerance = theSolid->GetTolerance();

  // Binary search to bring surface point to edge
  G4int maxItr = 100;
  G4int itr = 0;
  while ((high - low > tolerance || isIn != kSurface) && ++itr < maxItr) {
    mid = 0.5 * (low + high);
    localPos = originalPos + mid * vTan;

    // Modify localPos in place to surface if adjusting on curved surfaces
    if (curvedSurf == 1) AdjustToClosestSurfacePoint(localPos);
    isIn = theSolid->Inside(localPos);

    if (isIn == kSurface) low = mid; // Move out
    else high = mid; // Move in
  }

  if (verboseLevel>3) {
    G4cout << "G4CMPSolidUtils::AdjustToEdgePosition"
      << ": initialPos = " << originalPos
      << ", vTan = " << vTan
      << ", finalPos = " << localPos << G4endl;
  }
}

// Get edge position along vTan from localPos
G4ThreeVector G4CMPSolidUtils::GetEdgePosition(const G4ThreeVector& vTan,
                                               const G4ThreeVector& localPos,
                                               G4double high, const G4int curvedSurf) const {
  G4ThreeVector edgePosition = localPos;
  AdjustToEdgePosition(vTan, edgePosition, high, curvedSurf);
  return edgePosition;
}


// Reflect vTan against an edge
// Modifies vTan and surfNorm in place

void G4CMPSolidUtils::
ReflectAgainstEdge(G4ThreeVector& vTan, const G4ThreeVector& localPos,
                   G4ThreeVector& surfNorm) const {
  // Get normal of both surfaces at this edge point by stepping with normAdjust
  G4double normAdjust = 1*nm;
  G4double vTanMag = vTan.mag();
  G4ThreeVector norm1 = theSolid->SurfaceNormal(localPos);
  G4ThreeVector norm2 = theSolid->SurfaceNormal(localPos - normAdjust*norm1);
  // Try to fix norm1 if we didn't get the normals for a proper reflection
  norm1 = (norm1 * norm2 > 0) ? theSolid->SurfaceNormal(localPos - normAdjust*norm2) : norm1;

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

  if (verboseLevel>3) {
    G4cout << "G4CMPSolidUtils::ReflectAgainstEdge"
      << ": localPos = " << localPos
      << ", vTan_0 = " << vTan - 2*((vTan * -refNorm) * -refNorm)
      << ", edgeVector = " << edgeVec
      << ", refNorm = " << refNorm
      << ", norm1 = " << norm1
      << ", norm2 = " << norm2
      << ", vTan_f = " << vTan << G4endl;
  }
}
/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4PhononReflection.cc
/// \brief Implementation of the G4PhononReflection class
//
// This utility class handles general G4VSolid processes.
//
// $Id$
//
// 20250424  G4CMP-465 -- Create G4CMPSolidUtils class.

#include "G4CMPSolidUtils.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4VSolid.hh"


// Copy operations

G4CMPSolidUtils::G4CMPSolidUtils(const G4VSolid* solid)
  : theSolid(solid), verboseLevel(0) {;}

G4CMPSolidUtils&
G4CMPSolidUtils::operator=(const G4CMPSolidUtils& right) {
  theSolid = right.theSolid;
  verboseLevel = right.verboseLevel;
  return *this;
}
  
// Initialize the solid
  
void G4CMPSolidUtils::Initialize(const G4VSolid* solid) {
  theSolid = solid;
}


// Find the direction for minimum distance to surface
// Adjust theta0 and phi0 in place with best match
void G4CMPSolidUtils::
OptimizeSurfaceAdjustAngle(const G4ThreeVector& stepLocalPos, G4double& theta0,
                           G4double& phi0, const G4int angOption,
                           const G4double minDist) const {
  // Constants used for Golden Section searches
  G4double const tolerance = 1e-12*rad;
  G4double const gRatio = (std::sqrt(5) + 1) / 2;
  EInside const isIn = theSolid->Inside(stepLocalPos);

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
  G4double dist1 = (isIn == kInside) ? theSolid->DistanceToOut(stepLocalPos + dir1)
                                     : theSolid->DistanceToIn(stepLocalPos + dir1);
  G4double dist2 = (isIn == kInside) ? theSolid->DistanceToOut(stepLocalPos + dir2)
                                     : theSolid->DistanceToIn(stepLocalPos + dir2);

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
      dist1 = (isIn == kInside) ? theSolid->DistanceToOut(stepLocalPos + dir1)
                                : theSolid->DistanceToIn(stepLocalPos + dir1);
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
      dist2 = (isIn == kInside) ? theSolid->DistanceToOut(stepLocalPos + dir2)
                                : theSolid->DistanceToIn(stepLocalPos + dir2);
    }
  }
  G4double optimalAng = (x1 + x2) / 2;

  // Adjust angle in place
  if (angOption == 0) { theta0 = optimalAng; }
  else { phi0 = optimalAng; }
}


// Adjust to surface position closest to stepLocalPos without SurfaceNormal

void G4CMPSolidUtils::
AdjustToClosestSurfacePoint(G4ThreeVector& stepLocalPos) const {
  // Determine where you are in respect to solid
  EInside isIn = theSolid->Inside(stepLocalPos);
  if (isIn == kSurface) return;

  // Angles to be adjusted in place by OptimizeSurfaceAdjustAngle
  // Start at theta = pi/2 so "fit" can be determined by phi
  G4double bestTheta = pi / 2;
  G4double bestPhi = 0;

  G4double minDist = (isIn == kInside) ? theSolid->DistanceToOut(stepLocalPos)
                                       : theSolid->DistanceToIn(stepLocalPos);

  // Try to optimize phi first
  OptimizeSurfaceAdjustAngle(stepLocalPos, bestTheta, bestPhi, 1, minDist);
  OptimizeSurfaceAdjustAngle(stepLocalPos, bestTheta, bestPhi, 0, minDist);

  G4ThreeVector optDir(minDist*sin(bestTheta)*cos(bestPhi),
                       minDist*sin(bestTheta)*sin(bestPhi),
                       minDist*cos(bestTheta));

  // Only return valid positions on surface
  if (theSolid->Inside(stepLocalPos + optDir) == kSurface) {
    stepLocalPos += optDir;
  } else {
    stepLocalPos.set(kInfinity,kInfinity,kInfinity);
  }
}


// Do a binary search to find the closest point toward the edge along kTan
// Modifies stepLocalPos in place

void G4CMPSolidUtils::
AdjustToEdgePosition(const G4ThreeVector& kTan,
                     G4ThreeVector& stepLocalPos, G4double high,
                     const G4int curvedSurf) const {
  EInside isIn = theSolid->Inside(stepLocalPos);
  G4double low = 0.0*um;
  G4double mid = 0;
  G4ThreeVector originalPos = stepLocalPos;
  G4double tolerance = theSolid->GetTolerance();

  // Binary search to bring surface point to edge
  G4int maxItr = 100;
  G4int itr = 0;
  while ((high - low > tolerance || isIn != kSurface) && ++itr < maxItr) {
    mid = 0.5 * (low + high);
    stepLocalPos = originalPos + mid * kTan;

    // Modify stepLocalPos in place to surface if adjusting on curved surfaces
    if (curvedSurf == 1) AdjustToClosestSurfacePoint(stepLocalPos);
    isIn = theSolid->Inside(stepLocalPos);

    if (isIn == kSurface) low = mid; // Move out
    else high = mid; // Move in
  }

  if (verboseLevel>3) {
    G4double remDist = (theSolid->DistanceToIn(stepLocalPos) >=
                        theSolid->DistanceToOut(stepLocalPos))
                      ? theSolid->DistanceToIn(stepLocalPos)
                      : theSolid->DistanceToOut(stepLocalPos);
    G4cout << "G4CMPSolidUtils::AdjustToEdgePosition"
      << ": initialPos = " << originalPos
      << ", kTan = " << kTan
      << ", finalPos = " << stepLocalPos
      << ", remainingDist = " << remDist << G4endl;
  }
}


// Reflect kTan against an edge
// Modifies kTan and newNorm in place

void G4CMPSolidUtils::
ReflectAgainstEdge(G4ThreeVector& kTan,
                   const G4ThreeVector& stepLocalPos, G4ThreeVector& newNorm) const {
  // Get normal of both surfaces at this edge point by stepping with normAdjust
  G4double normAdjust = 1*nm;
  G4double kTanMag = kTan.mag();
  G4ThreeVector norm1 = theSolid->SurfaceNormal(stepLocalPos);
  G4ThreeVector norm2 = theSolid->SurfaceNormal(stepLocalPos - normAdjust*norm1);
  // Try to fix norm1 if we didn't get the normals for a proper reflection
  norm1 = (norm1 * norm2 > 0) ? theSolid->SurfaceNormal(stepLocalPos - normAdjust*norm2) : norm1;

  G4ThreeVector edgeVec(0,0,0);
  G4ThreeVector refNorm(0,0,0);

  // Only do reflection if at an edge
  if (norm1 * norm2 <= 0) {
    newNorm = (newNorm*norm1 > newNorm*norm2) ? norm1 : norm2;
    // Project kTan to be orthogonal to one normal
    (kTan -= newNorm * (kTan * newNorm)).setMag(kTanMag);

    // Get the edge vector
    edgeVec = norm1.cross(norm2).unit();

    // Find the normal for the reflection "surface"
    refNorm = (edgeVec.cross(newNorm)).unit();
    if (refNorm * kTan < 0) refNorm *= -1;

    // Reflect kTan against reflection "surface"
    kTan -= 2*((kTan * refNorm) * refNorm);
  }

  if (verboseLevel>3) {
    G4cout << "G4CMPSolidUtils::ReflectAgainstEdge"
      << ": stepLocalPos = " << stepLocalPos
      << ", kTan_0 = " << kTan - 2*((kTan * -refNorm) * -refNorm)
      << ", edgeVector = " << edgeVec
      << ", refNorm = " << refNorm
      << ", norm1 = " << norm1
      << ", norm2 = " << norm2
      << ", kTan_f = " << kTan << G4endl;
  }
}
/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPSolidUtils.hh
/// \brief Definition of the G4CMPSolidUtils class
//
// $Id$
//
// 20250424 G4CMP-465 -- Create G4CMPSolidUtils class.

#ifndef G4CMPSolidUtils_hh
#define G4CMPSolidUtils_hh 1

#include "globals.hh"
#include "G4ThreeVector.hh"

class G4VSolid;

class G4CMPSolidUtils {
  public:

  // Default constructor
  G4CMPSolidUtils() : theSolid(0), verboseLevel(0) {;}

  void Initialize(const G4VSolid* solid);

  void SetVerboseLevel(G4int verbose=0) { verboseLevel = verbose; }

  const G4VSolid* GetSolid() { return theSolid; }

  // Efficiently find direction with min distance to surface
  void OptimizeSurfaceAdjustAngle(const G4ThreeVector& stepLocalPos,
                                  G4double& theta0, G4double& phi0,
                                  const G4int angOption,
                                  const G4double minDist) const;

  // Modifies stepLocalPos in place
  void AdjustToClosestSurfacePoint(G4ThreeVector& stepLocalPos) const;

  void AdjustToEdgePosition(const G4ThreeVector& kTan,
                            G4ThreeVector& stepLocalPos,
                            G4double high, const G4int curvedSurf) const;

  // Modifies kTan and newNorm in place
  void ReflectAgainstEdge(G4ThreeVector& kTan,
                          const G4ThreeVector& stepLocalPos,
                          G4ThreeVector& newNorm) const;
  protected:
    G4CMPSolidUtils(const G4VSolid* solid);
    G4CMPSolidUtils& operator=(const G4CMPSolidUtils& right);

  private:
    const G4VSolid* theSolid;
    G4int verboseLevel;
};

#endif	/* G4CMPSolidUtils_hh */
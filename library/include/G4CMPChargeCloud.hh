/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPChargeCloud.hh
/// \brief Definition of the G4CMPChargeCloud class
///   Generates a collection of points distributed in a sphere around a
///   given center.  The distribution is uniform in angle, falls off
///   linearly with radius.  Size of cloud determined by lattice structure
///   (eight per unit cell for diamond).  If enclosing volume is specified
///   sphere will be "folded" inward at bounding surfaces.  If a touchable
///   is specified, both the enclosing volume and local-global transforms
///   will be extracted, and the final distribution returned in global
///   coordinates.
///
// $Id$

#ifndef G4CMPChargeCloud_hh
#define G4CMPChargeCloud_hh 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include <vector>

class G4LatticeLogical;
class G4LatticePhysical;
class G4VPhysicalVolume;
class G4VSolid;
class G4VTouchable;


class G4CMPChargeCloud {
public:
  G4CMPChargeCloud(const G4LatticeLogical* lat=0, const G4VSolid* solid=0);
  G4CMPChargeCloud(const G4LatticePhysical* lat, const G4VSolid* solid);
  explicit G4CMPChargeCloud(const G4VSolid* solid);

  G4CMPChargeCloud(const G4LatticeLogical* lat, const G4VPhysicalVolume* vol);
  G4CMPChargeCloud(const G4LatticePhysical* lat, const G4VPhysicalVolume* vol);
  explicit G4CMPChargeCloud(const G4VPhysicalVolume* vol);

  G4CMPChargeCloud(const G4LatticeLogical* lat, const G4VTouchable* touch);
  G4CMPChargeCloud(const G4LatticePhysical* lat, const G4VTouchable* touch);
  explicit G4CMPChargeCloud(const G4VTouchable* touch);

  virtual ~G4CMPChargeCloud();

  // Configure for operation
  void SetVerboseLevel(G4int vb) { verboseLevel=vb; }
  G4int GetVerboseLevel() const { return verboseLevel; }

  void SetLattice(const G4LatticeLogical* lat) { theLattice = lat; }
  void SetLattice(const G4LatticePhysical* lat);

  void SetTouchable(const G4VTouchable* touch);
  void UseVolume(const G4VPhysicalVolume* vol);
  void SetShape(const G4VSolid* solid) { theSolid = solid; }

  // Fill list of positions around specified center, within optional volume
  const std::vector<G4ThreeVector>& Generate(G4int npos, G4ThreeVector center);
  // NOTE:  Pass-by-value to allow global-to-local transform if necessary

  // Get previously filled list of positions
  const std::vector<G4ThreeVector>& GetCloud() const { return theCloud; }

  // Compute approximate radius of sphere to contain points
  G4double GetRadius(G4int npos) const;

  // Generate point randomly in sphere of given radius
  G4ThreeVector GeneratePoint(G4double rmax) const;

  // Adjust specified point to be inside volume
  void AdjustToVolume(G4ThreeVector& point) const;

protected:
  G4int verboseLevel;			// Diagnostic messages
  const G4LatticeLogical* theLattice;	// For crystal structure
  const G4VSolid* theSolid;		// For bounding surfaces
  const G4VTouchable* theTouchable;	// For local-global coordinates

private:
  std::vector<G4ThreeVector> theCloud;	// Buffer to carry generated points
};

#endif	/* G4CMPChargeCloud_hh */

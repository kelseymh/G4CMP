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
//
// 20170925  Add direct access to individual positions in cloud, binning
// 20180831  Fix compiler warning on GetPositionBin()

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
  explicit G4CMPChargeCloud(const G4VPhysicalVolume* vol);
  explicit G4CMPChargeCloud(const G4VTouchable* touch);

  virtual ~G4CMPChargeCloud();

  // Configure for operation
  void SetVerboseLevel(G4int vb) { verboseLevel=vb; }
  G4int GetVerboseLevel() const { return verboseLevel; }

  void SetLattice(const G4LatticeLogical* lat);
  void SetLattice(const G4LatticePhysical* lat);
  const G4LatticeLogical* GetLattice() const { return theLattice; }

  void SetTouchable(const G4VTouchable* touch);
  const G4VTouchable* GetTouchable() const { return theTouchable; }

  void UseVolume(const G4VPhysicalVolume* vol);
  void SetShape(const G4VSolid* solid) { theSolid = solid; }
  const G4VSolid* GetShape() const { return theSolid; }

  // Fill list of positions around specified center, within optional volume
  virtual const std::vector<G4ThreeVector>& 
  Generate(G4int npos, const G4ThreeVector& center);

  // Get constituents and parameters of generated distribution
  const std::vector<G4ThreeVector>& GetCloud() const { return theCloud; }
  const G4ThreeVector& GetPosition(G4int i) const { return theCloud[i]; }

  const std::vector<G4int>& GetCloudBins() const { return theCloudBins; }
  G4int GetPositionBin(G4int i) const { return theCloudBins[i]; }
  G4ThreeVector GetBinCenter(G4int ibin) const;

  G4double GetRadius() const { return cloudRadius; }
  const G4ThreeVector& GetCenter() const { return localCenter; }

  // Compute maximum radius of sphere to contain points
  virtual G4double ComputeRadius(G4int npos) const;

  // Generate point randomly in sphere of given radius
  virtual G4ThreeVector GeneratePoint(G4double rmax) const;

  // Adjust specified point to be inside volume
  void AdjustToVolume(G4ThreeVector& point) const;

protected:
  G4int verboseLevel;			// Diagnostic messages
  const G4LatticeLogical* theLattice;	// For crystal structure
  const G4VSolid* theSolid;		// For bounding surfaces
  const G4VTouchable* theTouchable;	// For local-global coordinates
  G4double avgLatticeSpacing;		// Cubic approximation from lattice
  G4double radiusScale;			// Cloud radius per e/h pair (cbrt)
  G4double binSpacing;			// Bin size for primary clumping

  // Convert local position to bin index (pass-by-value for use as temporary)
  G4int GetBinIndex(G4ThreeVector localPos) const;

private:
  std::vector<G4ThreeVector> theCloud;	// Buffer to carry generated points
  G4double cloudRadius;			// Radius used to generate distribution
  G4ThreeVector localCenter;		// Local center point of distribution
  std::vector<G4int> theCloudBins;	// Buffer for bin indices at points
};

#endif	/* G4CMPChargeCloud_hh */

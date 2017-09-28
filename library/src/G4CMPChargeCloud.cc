/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPChargeCloud.cc
/// \brief Implementation of the G4CMPChargeCloud class
///   Generates a collection of points distributed in a sphere around a
///   given center.  The distribution is uniform in angle, falls off
///   linearly with radius.  Size of cloud determined by lattice structure
///   (eight per unit cell for diamond).  If enclosing volume is specified
///   sphere will be "folded" inward at bounding surfaces.
///
// $Id$

#include "G4CMPChargeCloud.hh"
#include "G4CMPGeometryUtils.hh"
#include "G4LatticeLogical.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalVolume.hh"
#include "G4RandomDirection.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"
#include "G4VTouchable.hh"
#include "Randomize.hh"
#include <math.h>


// Constructors and destructor

G4CMPChargeCloud::G4CMPChargeCloud(const G4LatticeLogical* lat,
				   const G4VSolid* solid)
  : verboseLevel(0), theLattice(0), theSolid(solid), theTouchable(nullptr),
    avgLatticeSpacing(0.), radiusScale(0.), binSpacing(0.), cloudRadius(0.) {
  SetLattice(lat);
}

G4CMPChargeCloud::G4CMPChargeCloud(const G4LatticePhysical* lat,
				   const G4VSolid* solid)
  : G4CMPChargeCloud(lat->GetLattice(), solid) {;}

G4CMPChargeCloud::G4CMPChargeCloud(const G4VSolid* solid)
  : G4CMPChargeCloud((const G4LatticeLogical*)nullptr, solid) {;}

G4CMPChargeCloud::G4CMPChargeCloud(const G4VPhysicalVolume* vol)
  : G4CMPChargeCloud() {
  UseVolume(vol);
}

G4CMPChargeCloud::G4CMPChargeCloud(const G4VTouchable* touch)
  : G4CMPChargeCloud() {
  SetTouchable(touch);
}

G4CMPChargeCloud::~G4CMPChargeCloud() {;}


// Store lattice and compute derived quantities

void G4CMPChargeCloud::SetLattice(const G4LatticePhysical* lat) {
  SetLattice(lat ? lat->GetLattice() : nullptr);
}

void G4CMPChargeCloud::SetLattice(const G4LatticeLogical* lat) {
  theLattice = lat;

  // Equivalent cubic unit cell dimensions regardless of crystal group
  avgLatticeSpacing = !lat ? 1.*nm : cbrt(lat->GetBasis(0).mag() *
					  lat->GetBasis(1).mag() *
					  lat->GetBasis(2).mag());

  radiusScale = avgLatticeSpacing/4.;	// 8 atoms per unit cell, r = a0/2.

  binSpacing = 4.*avgLatticeSpacing;	// Box sizes for collecting primaries
}


// Store touchable (transform) and volume shape

void G4CMPChargeCloud::SetTouchable(const G4VTouchable* touch) {
  UseVolume(touch ? touch->GetVolume() : nullptr);
  theTouchable = touch;
}


// Extract shape from placement volume

void G4CMPChargeCloud::UseVolume(const G4VPhysicalVolume* vol) {
  theTouchable = nullptr;
  SetShape(vol ? vol->GetLogicalVolume()->GetSolid() : nullptr);

  SetLattice(G4LatticeManager::Instance()->GetLattice(vol));
}


// Fill list of positions around specified center, within optional volume
// If user specified G4VTouchable, coordinates are all global

const std::vector<G4ThreeVector>& 
G4CMPChargeCloud::Generate(G4int npos, const G4ThreeVector& center) {
  localCenter = center;
  if (theTouchable) G4CMP::RotateToLocalPosition(theTouchable, localCenter);

  cloudRadius = ComputeRadius(npos);	// Radius of cloud for average density

  if (verboseLevel) {
    G4cout << "G4CMPChargeCloud::Generate " << npos << " @ " << center
	   << " radius " << cloudRadius/nm << " bins "
	   << binSpacing/nm << " nm" << G4endl;
  }

  theCloud.clear();
  theCloud.reserve(npos);

  theCloudBins.clear();
  theCloudBins.reserve(npos);

  for (G4int i=0; i<npos; i++) {
    theCloud.push_back(GeneratePoint(cloudRadius)+localCenter);
    if (theSolid) AdjustToVolume(theCloud.back());	// Checkout boundaries

    theCloudBins.push_back(GetBinIndex(theCloud.back()));

    if (theTouchable)
      G4CMP::RotateToGlobalPosition(theTouchable, theCloud.back());

    if (verboseLevel>2) {
      G4cout << " point " << i << " @ " << theCloud.back() << " in bin "
	     << theCloudBins.back() << G4endl;
    }
  }

  return theCloud;
}


// Compute radius of cloud for average density matching unit cell

G4double G4CMPChargeCloud::ComputeRadius(G4int npos) const {
  // Charge cloud must be at least a full unit cell in size
  G4double rmax = std::max(avgLatticeSpacing, radiusScale*cbrt(npos));

  if (verboseLevel>1) {
    G4cout << "G4CMPChargeCloud unit cell " << avgLatticeSpacing/nm << " nm "
	   << " radius " << rmax/nm << " nm for " << npos << " charges"
	   << G4endl;
  }

  return rmax;
}


// Generate random point in specified sphere

G4ThreeVector G4CMPChargeCloud::GeneratePoint(G4double rmax) const {
  G4double rndm = G4UniformRand();
  G4double r = rmax*(1.-sqrt(1.-rndm*rndm));	// Linear from r=0 to rmax

  return r*G4RandomDirection();
}


// Adjust specified point to be inside volume

void G4CMPChargeCloud::AdjustToVolume(G4ThreeVector& point) const {
  if (!theSolid) return;			// Avoid unnecessary work

  // Keep adjusting point until interior to volume
  G4ThreeVector norm;
  G4double dist = 0.;

  G4int ntries = 100;				// Avoid infinite loops
  while (--ntries > 0 && theSolid->Inside(point) == kOutside) {
    norm = theSolid->SurfaceNormal(point);
    dist = theSolid->DistanceToIn(point);

    point -= 2.*dist*norm;	// Reflect point through boundary surface
  }
}


// Conversions between local position and bin index (ijk)
// Binning is done in units of 4 lattice spacings, centered on (0,0,0)

G4int G4CMPChargeCloud::GetBinIndex(G4ThreeVector pos) const {
  // NOTE:  Argument passed by value for use in computations below

  // Convert local position in volume to bin index (w/(0,0,0) at bin center)
  pos -= localCenter;
  pos += G4ThreeVector(cloudRadius+binSpacing/2., cloudRadius+binSpacing/2.,
		       cloudRadius+binSpacing/2.);

  pos /= binSpacing;

  // Bin index (ijk) allows up to 64 billion unit cells in cloud
  return (int(pos.z())*1000 + int(pos.y()))*1000 + int(pos.x());
}

G4ThreeVector G4CMPChargeCloud::GetBinCenter(G4int ibin) const {
  G4ThreeVector pbin((ibin%1000 - 0.5)*binSpacing - cloudRadius,
		     ((ibin/1000)%1000 - 0.5)*binSpacing - cloudRadius,
		     (ibin/1000000 - 0.5)*binSpacing - cloudRadius);
  pbin += localCenter;

  return pbin;
}

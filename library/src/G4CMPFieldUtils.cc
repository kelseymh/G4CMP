/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
// File: G4CMPFieldUtils.cc
//
// Description: Free standing helper functions for electric field access
//
// 20180622  Michael Kelsey

#include "G4CMPFieldUtils.hh"
#include "G4CMPGeometryUtils.hh"
#include "G4CMPLocalElectroMagField.hh"
#include "G4CMPMeshElectricField.hh"
#include "G4ElectroMagneticField.hh"
#include "G4Field.hh"
#include "G4FieldManager.hh"
#include "G4LogicalVolume.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4UniformElectricField.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"
#include "G4VTouchable.hh"


// Get electric field at specified position (track or step)

namespace {
  G4ThreeVector origin(0.,0.,0.);	// For convenience below
}

G4ThreeVector G4CMP::GetFieldAtPosition(const G4Step& step) {
  return GetFieldAtPosition(*(step.GetTrack()));
}

G4ThreeVector G4CMP::GetFieldAtPosition(const G4Track& track) {
  return GetFieldAtPosition(track.GetTouchable(), track.GetPosition());
}


// Get field using _global_ coordinates, for specified volume

G4ThreeVector G4CMP::GetFieldAtPosition(const G4VTouchable* touch,
					const G4ThreeVector& pos) {
  const G4LogicalVolume* vol = (touch ? touch->GetVolume()->GetLogicalVolume()
				: GetVolumeAtPoint(pos)->GetLogicalVolume());
  const G4FieldManager* fMan = vol->GetFieldManager();

  if (!fMan || !fMan->DoesFieldExist()) return origin;

  G4double position[4] = { pos[0], pos[1], pos[2], 0. };
  G4double fieldVal[6];
  fMan->GetDetectorField()->GetFieldValue(position, fieldVal);

  return G4ThreeVector(fieldVal[3], fieldVal[4], fieldVal[5]);
}


// Get potential (uniform or mesh field) at _global_ coordinate

G4double G4CMP::GetPotentialAtPosition(const G4Step& step) {
  return GetPotentialAtPosition(*(step.GetTrack()));
}

G4double G4CMP::GetPotentialAtPosition(const G4Track& track) {
  return GetPotentialAtPosition(track.GetTouchable(), track.GetPosition());
}

// Get potential at starting position of track

G4double G4CMP::GetPotentialAtVertex(const G4Track& track) {
  return GetPotentialAtPosition(track.GetTouchable(),
				track.GetVertexPosition());
}


// Get potential using _global_ coordinates, for specified volume

G4double G4CMP::GetPotentialAtPosition(const G4VTouchable* touch,
				       const G4ThreeVector& pos) {
  const G4LogicalVolume* vol = (touch ? touch->GetVolume()->GetLogicalVolume()
				: GetVolumeAtPoint(pos)->GetLogicalVolume());

  G4double position[4] = { pos[0], pos[1], pos[2], 0. };

  // G4CMP special fields can return potential directly
  if (GetLocalField(vol)) return GetLocalField(vol)->GetPotential(position);
  if (GetMeshField(vol)) return GetMeshField(vol)->GetPotential(position);

  // Uniform field can be "integrated" across volume, assume V=0 at midplane
  const G4UniformElectricField* ufield = GetUniformField(vol);
  if (ufield) {
    G4double be[6];                             // Bx, By, Bz, Ex, Ey, Ez
    ufield->GetFieldValue(position, be);
    G4ThreeVector evec(be[3],be[4],be[5]);      // Field is same everywhere

    G4VSolid* shape = vol->GetSolid();
    G4ThreeVector e0 = evec.unit();
    G4double toVpos = shape->DistanceToOut(pos, -e0);
    G4double toVneg = shape->DistanceToOut(pos, e0);

    return 0.5*(toVpos-toVneg)*evec.mag();	// [-V/2,V/2] interpolation
  }

  // Arbitrary field configurations must be integrated
  return 0.;
}


// Estimate bias across volume using input position as reference
// NOTE:  This doesn't integrate field, just uses local magnitude, direction

G4double G4CMP::GetBiasThroughPosition(const G4VTouchable* touch,
				       const G4ThreeVector& pos) {
  G4ThreeVector field = GetFieldAtPosition(touch, pos);
  if (field.isNear(origin)) return 0.;		// No field, no bias estimate

  const G4LogicalVolume* vol = (touch ? touch->GetVolume()->GetLogicalVolume()
				: GetVolumeAtPoint(pos)->GetLogicalVolume());
  G4VSolid* shape = vol->GetSolid();

  // Thickness of volume through point along local field direction
  // NOTE: This is only valid for uniform or near-uniform fields
  G4ThreeVector e0 = field.unit();
  G4double thick = shape->DistanceToOut(pos,-e0) + shape->DistanceToOut(pos,e0);

  // Arbitrary field configurations should be integrated along lines
  return field.mag() * thick;
}


// Return different field handlers if available for volume

const G4CMPLocalElectroMagField* 
G4CMP::GetLocalField(const G4LogicalVolume* vol) {
  const G4FieldManager* fMan = vol->GetFieldManager();
  if (!fMan || !fMan->DoesFieldExist()) return 0;

  const G4Field* baseField = fMan->GetDetectorField();
  return dynamic_cast<const G4CMPLocalElectroMagField*>(baseField);
}
  
const G4CMPMeshElectricField* 
G4CMP::GetMeshField(const G4LogicalVolume* vol) {
  const G4Field* field = 0;

  if (GetLocalField(vol)) field = GetLocalField(vol)->GetLocalField();
  else {
    const G4FieldManager* fMan = vol->GetFieldManager();
    if (fMan && fMan->DoesFieldExist()) field = fMan->GetDetectorField();
  }

  return dynamic_cast<const G4CMPMeshElectricField*>(field);
}

const G4UniformElectricField* 
G4CMP::GetUniformField(const G4LogicalVolume* vol) {
  const G4Field* field = 0;

  if (GetLocalField(vol)) field = GetLocalField(vol)->GetLocalField();
  else {
    const G4FieldManager* fMan = vol->GetFieldManager();
    if (fMan && fMan->DoesFieldExist()) field = fMan->GetDetectorField();
  }

  return dynamic_cast<const G4UniformElectricField*>(field);
}

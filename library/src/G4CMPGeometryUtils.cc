/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
// File: G4CMPGeometryUtils.cc
//
// Description: Free standing helper functions for geometry based calculations.
//
// 20161107  Rob Agnese

#include "G4CMPGeometryUtils.hh"

#include "G4CMPGlobalLocalTransformStore.hh"

G4ThreeVector G4CMP::GetLocalDirection(const G4VPhysicalVolume* pv,
                                       const G4ThreeVector& dir) {
  return G4CMPGlobalLocalTransformStore::ToLocal(pv).TransformAxis(dir);
}

G4ThreeVector G4CMP::GetLocalPosition(const G4VPhysicalVolume* pv,
                                      const G4ThreeVector& pos) {
  return G4CMPGlobalLocalTransformStore::ToLocal(pv).TransformPoint(pos);
}

G4ThreeVector G4CMP::GetGlobalDirection(const G4VPhysicalVolume* pv,
                                        const G4ThreeVector& dir) {
  return G4CMPGlobalLocalTransformStore::ToGlobal(pv).TransformAxis(dir);
}

G4ThreeVector G4CMP::GetGlobalPosition(const G4VPhysicalVolume* pv,
                                       const G4ThreeVector& pos) {
  return G4CMPGlobalLocalTransformStore::ToGlobal(pv).TransformPoint(pos);
}

void G4CMP::RotateToLocalDirection(const G4VPhysicalVolume* pv,
                                   G4ThreeVector& dir) {
  G4CMPGlobalLocalTransformStore::ToLocal(pv).ApplyAxisTransform(dir);
}

void G4CMP::RotateToLocalPosition(const G4VPhysicalVolume* pv,
                                  G4ThreeVector& pos) {
  G4CMPGlobalLocalTransformStore::ToLocal(pv).ApplyPointTransform(pos);
}

void G4CMP::RotateToGlobalDirection(const G4VPhysicalVolume* pv,
                                    G4ThreeVector& dir) {
  G4CMPGlobalLocalTransformStore::ToGlobal(pv).ApplyAxisTransform(dir);
}

void G4CMP::RotateToGlobalPosition(const G4VPhysicalVolume* pv,
                                   G4ThreeVector& pos) {
  G4CMPGlobalLocalTransformStore::ToGlobal(pv).ApplyPointTransform(pos);
}

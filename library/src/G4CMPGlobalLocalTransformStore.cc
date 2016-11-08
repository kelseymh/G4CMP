/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
// File: G4CMPGlobalLocalTransformStore.cc
//
// Description: Singleton cache for G4AffineTransforms that are frequently
//      used for coordinate transformations between the global coordinate
//      system (used by internal G4 tracking) and the local-to-solid coordinates
//      (used by MeshElectricField and most lattice-aware physics processes).
//
// 20161102  Rob Agnese

#include "G4CMPGlobalLocalTransformStore.hh"

#include "G4VPhysicalVolume.hh"

const G4AffineTransform&
G4CMPGlobalLocalTransformStore::ToLocal(const G4VPhysicalVolume* pv) {
  return Instance().GetOrBuildTransforms(pv).globalToLocal;
}

const G4AffineTransform&
G4CMPGlobalLocalTransformStore::ToGlobal(const G4VPhysicalVolume* pv) {
  return Instance().GetOrBuildTransforms(pv).localToGlobal;
}

void
G4CMPGlobalLocalTransformStore::Reset() {
  Instance().cache.clear();
}

const G4CMPGlobalLocalTransformStore::Transforms&
G4CMPGlobalLocalTransformStore::GetOrBuildTransforms(const G4VPhysicalVolume* pv) {
  if (!pv) {
    G4Exception("G4CMPGlobalLocalTransformStore::GetOrBuildTransforms",
                "trans001", FatalErrorInArgument,
                "PhysicalVolume pointer is null.");
  }

  if (Instance().cache.count(pv) == 0) {
    auto localToGlobal = G4AffineTransform(pv->GetRotation(),
                                           pv->GetTranslation());
    return Instance().cache[pv] = Transforms { localToGlobal,
                                               localToGlobal.Inverse() };
  }

  return Instance().cache[pv];
}

G4CMPGlobalLocalTransformStore&
G4CMPGlobalLocalTransformStore::Instance() {
  static G4CMPGlobalLocalTransformStore instance;
  return instance;
}

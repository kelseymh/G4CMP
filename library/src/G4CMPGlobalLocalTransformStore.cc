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

#include "G4AutoLock.hh"
namespace {
  G4Mutex mut = G4MUTEX_INITIALIZER;
}

G4AffineTransform&
G4CMPGlobalLocalTransformStore::ToLocal(const G4VPhysicalVolume* pv) {
  return Instance().GetOrBuildTransforms(pv).globalToLocal;
}

G4AffineTransform&
G4CMPGlobalLocalTransformStore::ToGlobal(const G4VPhysicalVolume* pv) {
  return Instance().GetOrBuildTransforms(pv).localToGlobal;
}

void
G4CMPGlobalLocalTransformStore::Reset() {
  Instance().cache.clear();
}

G4CMPGlobalLocalTransformStore::Transforms&
G4CMPGlobalLocalTransformStore::GetOrBuildTransforms(const G4VPhysicalVolume* pv) {
  //if (!pv) return Transforms(); // Identity for both

  if (Instance().cache.count(pv) == 0) {
    auto localToGlobal = G4AffineTransform(pv->GetRotation(),
                                           pv->GetTranslation());
    G4AutoLock lock(&mut);
    return Instance().cache[pv] = Transforms { localToGlobal,
                                               localToGlobal.Inverse() };
    // unlock
  }

  return Instance().cache[pv];
}

G4CMPGlobalLocalTransformStore&
G4CMPGlobalLocalTransformStore::Instance() {
  static G4CMPGlobalLocalTransformStore instance;
  return instance;
}

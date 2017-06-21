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
// 20170605  Pass touchable from track, not just local TOUCH

#include "G4CMPGlobalLocalTransformStore.hh"

#include "G4NavigationHistory.hh"
#include "G4VTouchable.hh"


G4CMPGlobalLocalTransformStore&
G4CMPGlobalLocalTransformStore::Instance() {
  static G4CMPGlobalLocalTransformStore instance;
  return instance;
}

const G4AffineTransform&
G4CMPGlobalLocalTransformStore::ToLocal(const G4VTouchable* touch) {
  return Instance().GetOrBuildTransforms(touch).globalToLocal;
}

const G4AffineTransform&
G4CMPGlobalLocalTransformStore::ToGlobal(const G4VTouchable* touch) {
  return Instance().GetOrBuildTransforms(touch).localToGlobal;
}

void
G4CMPGlobalLocalTransformStore::Reset() {
  Instance().cache.clear();
}

const G4CMPGlobalLocalTransformStore::Transforms&
G4CMPGlobalLocalTransformStore::GetOrBuildTransforms(const G4VTouchable* touch) {
  if (!touch) {
    G4Exception("G4CMPGlobalLocalTransformStore::GetOrBuildTransforms",
                "trans001", FatalErrorInArgument,
                "Touchable pointer is null.");
  }

  uintptr_t thash = Hash(touch);
  if (Instance().cache.count(thash) == 0) {
    const G4AffineTransform& lToG = touch->GetHistory()->GetTransform(0);
    Instance().cache[thash] = Transforms { lToG, lToG.Inverse() };
  }

  return Instance().cache[thash];
}

// Convert touchable volume chain to unique identifier
uintptr_t 
G4CMPGlobalLocalTransformStore::Hash(const G4VTouchable* touch) const {
  if (!touch) {
    G4Exception("G4CMPGlobalLocalTransformStore::GetOrBuildTransforms",
                "trans001", FatalErrorInArgument,
                "Touchable pointer is null.");
  }

  static const uintptr_t prime = 12764787846358441471UL;	// 2^64 / Phi
  uintptr_t h = 1;
  for (G4int d=0; d<touch->GetHistoryDepth(); d++) {
    // NOTE:  Mapping pointer address to integer must be done explicitly
    uintptr_t pv = (uintptr_t)(void*)(touch->GetVolume(d));
    h *= prime + 2*(pv%prime);
    h %= prime;
  }

  return h/2;
}

/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef G4CMPGlobalLocalTransformStore_hh
#define G4CMPGlobalLocalTransformStore_hh 1

// $Id$
// File: G4CMPGlobalLocalTransformStore.hh
//
// Description: Singleton cache for G4AffineTransforms that are frequently
//      used for coordinate transformations between the global coordinate
//      system (used by internal G4 tracking) and the local-to-solid coordinates
//      (used by MeshElectricField and most lattice-aware physics processes).
//
// 20161102  Rob Agnese
// 20170605  Pass touchable from track, not just local PV
// 20200519  Convert to thread-local singleton (for use by worker threads)

#include "G4AffineTransform.hh"
#include "G4ThreadLocalSingleton.hh"
#include <cstdint>
#include <unordered_map>

class G4VTouchable;


class G4CMPGlobalLocalTransformStore {
public:
  static const G4AffineTransform& ToLocal(const G4VTouchable*);
  static const G4AffineTransform& ToGlobal(const G4VTouchable*);
  
  static void Reset();
  
private:
  friend class G4ThreadLocalSingleton<G4CMPGlobalLocalTransformStore>;

  G4CMPGlobalLocalTransformStore() = default;
  static G4CMPGlobalLocalTransformStore& Instance();
  
  struct Transforms {
    G4AffineTransform localToGlobal;
    G4AffineTransform globalToLocal;
  };
  
  std::unordered_map<uintptr_t, Transforms> cache;
  uintptr_t Hash(const G4VTouchable*) const;
  const Transforms& GetOrBuildTransforms(const G4VTouchable*);
};

#endif

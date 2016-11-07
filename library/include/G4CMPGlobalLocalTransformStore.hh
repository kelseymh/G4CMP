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

#include "G4AffineTransform.hh"
#include "unordered_map"

class G4VPhysicalVolume;

class G4CMPGlobalLocalTransformStore {
  public:
    static G4AffineTransform& ToLocal(const G4VPhysicalVolume*);
    static G4AffineTransform& ToGlobal(const G4VPhysicalVolume*);

    static void Reset();

  private:
    G4CMPGlobalLocalTransformStore() = default;
    static G4CMPGlobalLocalTransformStore& Instance();

    struct Transforms {
      G4AffineTransform localToGlobal;
      G4AffineTransform globalToLocal;
    };

    std::unordered_map<const G4VPhysicalVolume*, Transforms> cache;
    Transforms& GetOrBuildTransforms(const G4VPhysicalVolume*);
};

#endif

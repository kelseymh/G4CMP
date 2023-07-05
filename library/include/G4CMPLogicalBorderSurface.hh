/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/
// $Id$
// G4CMPLogicalBorderSurface
//
// Class description:
//
// A Logical Surface class for surfaces defined by the boundary
// of two physical volumes.
//
// Adapted from G4LogicalBorderSurface for phonon/charge carrier transport

#ifndef G4CMPLogicalBorderSurface_h
#define G4CMPLogicalBorderSurface_h 1

#include "G4LogicalSurface.hh"
#include <map>
#include <utility>

class G4VPhysicalVolume;
class G4CMPLogicalBorderSurface;

using G4CMPLogicalBorderKey = std::pair<const G4VPhysicalVolume*,
					const G4VPhysicalVolume*>;
using G4CMPLogicalBorderTable = std::map<G4CMPLogicalBorderKey,G4CMPLogicalBorderSurface*>;

class G4CMPLogicalBorderSurface : public G4LogicalSurface {
public:
  G4CMPLogicalBorderSurface(const G4String& name,
			    G4VPhysicalVolume* vol1,
			    G4VPhysicalVolume* vol2,
			    G4SurfaceProperty* surfaceProperty);
  virtual ~G4CMPLogicalBorderSurface() {;}

  // Assignment and copying not allowed.
  G4CMPLogicalBorderSurface(const G4CMPLogicalBorderSurface&) = delete;
  G4CMPLogicalBorderSurface& operator=(const G4CMPLogicalBorderSurface&) = delete;

  // Comparison operators match by pointer, not by contents!
  G4bool operator==(const G4CMPLogicalBorderSurface &right) const {
    return (this == &right);
  }

  G4bool operator!=(const G4CMPLogicalBorderSurface &right) const {
    return (this != &right);
  }

  // Get inbound and outbound volumes at surface
  const G4VPhysicalVolume* GetVolume1() const { return Volume1; }
  const G4VPhysicalVolume* GetVolume2() const { return Volume2; }

  // Reassign object in lookup table
  void SetVolumes(G4VPhysicalVolume* vol1, G4VPhysicalVolume* vol2);

  // Access lookup table mapping surfaces to volumes
  static G4CMPLogicalBorderSurface* GetSurface(const G4VPhysicalVolume* vol1,
					       const G4VPhysicalVolume* vol2);

  static const G4CMPLogicalBorderTable* GetSurfaceTable();
  static void CleanSurfaceTable();
  static size_t GetNumberOfSurfaces() {
    return surfaceTable ? surfaceTable->size() : 0;
  }

  // Report contents of lookup table for diagnostics
  static void DumpInfo(); // const 

private:
  // Used by SetVolumes() for entering or replacing object in registry
  void RemoveFromTable();		// Call _before_ changing volumes
  void AddToTable();			// Call _after_ changing volumes

private:
  G4VPhysicalVolume* Volume1;      // Physical Volume on side 1 (inbound)
  G4VPhysicalVolume* Volume2;      // Physical Volume on side 2 (outbound)

  // Table of BorderSurfaces -- created only in master thread
  static G4CMPLogicalBorderTable *surfaceTable;
};

#endif /* G4CMPLogicalBorderSurface_h */


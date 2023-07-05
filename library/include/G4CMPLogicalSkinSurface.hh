/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/
// $Id$
// G4CMPLogicalSkinSurface
//
// Class description:
//
// A Logical Surface class for the surface surrounding a single logical
// volume.
//
// Adapted from G4LogicalSkinSurface for phonon/charge carrier transport

#ifndef G4CMPLogicalSkinSurface_h
#define G4CMPLogicalSkinSurface_h 1

#include "G4LogicalSurface.hh"
#include <map>

class G4LogicalVolume;
class G4CMPLogicalSkinSurface;

using G4CMPLogicalSkinTable = std::map<const G4LogicalVolume*,
				       G4CMPLogicalSkinSurface*>;

class G4CMPLogicalSkinSurface : public G4LogicalSurface {
public:
  G4CMPLogicalSkinSurface(const G4String& name, G4LogicalVolume* vol,
		       G4SurfaceProperty* surfaceProperty);
  virtual ~G4CMPLogicalSkinSurface() {;}

  // Assignment and copying not allowed.
  G4CMPLogicalSkinSurface(const G4CMPLogicalSkinSurface&) = delete;
  G4CMPLogicalSkinSurface& operator=(const G4CMPLogicalSkinSurface&) = delete;

  // Comparison operators match by pointer, not by contents!
  G4bool operator==(const G4CMPLogicalSkinSurface &right) const {
    return (this == &right);
  }

  G4bool operator!=(const G4CMPLogicalSkinSurface &right) const {
    return (this != &right);
  }

  // Adjust configuration of object and lookup table
  const G4LogicalVolume* GetLogicalVolume() const { return LogVolume; }
  void SetLogicalVolume(G4LogicalVolume* vol);

  // Access lookup table mapping surfaces to volumes
  static G4CMPLogicalSkinSurface* GetSurface(const G4LogicalVolume* vol);

  static const G4CMPLogicalSkinTable* GetSurfaceTable();
  static void CleanSurfaceTable();
  static size_t GetNumberOfSurfaces() {
    return surfaceTable ? surfaceTable->size() : 0;
  }

  // Report contents of lookup table for diagnostics
  static void DumpInfo();
  
private:
  // Used by SetVolumes() for entering or replacing object in registry
  void RemoveFromTable();		// Call _before_ changing volumes
  void AddToTable();			// Call _after_ changing volumes

private:
  G4LogicalVolume* LogVolume;      // Logical Volume enclosed by skin

  // Table of SkinSurfaces -- created only in master thread
  static G4CMPLogicalSkinTable *surfaceTable;
};

#endif /* G4CMPLogicalSkinSurface_h */


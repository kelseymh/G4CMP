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
//
// 20220331  G4CMP-294:  Protect against null surfaceTable pointer

#include"G4CMPLogicalSkinSurface.hh"
#include "G4LogicalVolume.hh"


// Shared singleton -- instantiated during geom building on master only

G4CMPLogicalSkinTable* G4CMPLogicalSkinSurface::surfaceTable = 0;

const G4CMPLogicalSkinTable* G4CMPLogicalSkinSurface::GetSurfaceTable() {
  if (!surfaceTable) surfaceTable = new G4CMPLogicalSkinTable;
  return surfaceTable;
}


// Constructor

G4CMPLogicalSkinSurface::
G4CMPLogicalSkinSurface(const G4String& name, G4LogicalVolume* vol,
			G4SurfaceProperty* surfaceProperty)
  : G4LogicalSurface(name, surfaceProperty), LogVolume(0) {
  if (!surfaceTable) surfaceTable = new G4CMPLogicalSkinTable;
  SetLogicalVolume(vol);		// Also registers into map
}


// Reassign object in lookup table

void G4CMPLogicalSkinSurface::SetLogicalVolume(G4LogicalVolume* vol) {
  if (!vol) {
    G4Exception("G4CMPLogicalSkinSurface", "G4CMPSurf001",
		FatalException, "Null pointer passed for logical volume.");
    return;
  }

  RemoveFromTable();
  LogVolume = vol;
  AddToTable();
}


// Used by SetVolumes() for entering or replacing object in registry

void G4CMPLogicalSkinSurface::RemoveFromTable() {
  if (surfaceTable) surfaceTable->erase(LogVolume);
}

void G4CMPLogicalSkinSurface::AddToTable() {
  if (!surfaceTable) GetSurfaceTable();
  (*surfaceTable)[LogVolume] = this;
}


// Access lookup table mapping surfaces to volumes

G4CMPLogicalSkinSurface* 
G4CMPLogicalSkinSurface::GetSurface(const G4LogicalVolume* vol) {
  if (!surfaceTable) return 0;

  G4CMPLogicalSkinTable::iterator entry = surfaceTable->find(vol);
  return (entry != surfaceTable->end() ? entry->second : 0);
}


void G4CMPLogicalSkinSurface::CleanSurfaceTable() {
  if (!surfaceTable) return;		// No table defined, nothing to clear

  // Delete all defined surfaces, then remove table entries
  for (auto pos=surfaceTable->cbegin(); pos!=surfaceTable->cend(); ++pos) {
    if (pos->second) delete pos->second;
  }

  surfaceTable->clear();
}


// Report contents of registry for diagnostics

void G4CMPLogicalSkinSurface::DumpInfo() {
  G4cout << "***** G4CMPLogicalSkinSurface Table : Nb of Surfaces = "
         << GetNumberOfSurfaces() << " *****" << G4endl;

  if (!surfaceTable) return;

  for (auto pos=surfaceTable->cbegin(); pos!=surfaceTable->cend(); ++pos) {
    G4cout << pos->second->GetName() << " : " << G4endl
	   << " Surface of volume " << pos->second->GetLogicalVolume()->GetName()
	   << G4endl;
  }

  G4cout << G4endl;
}

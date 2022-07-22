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
//
// 20220331  G4CMP-294:  Protect against null surfaceTable pointer
// 20220718  G4CMP-306:  Ensure build compatibility with Geant4-10.04

#include "G4CMPLogicalBorderSurface.hh"
#include "G4ExceptionSeverity.hh"
#include "G4LogicalSurface.hh"
#include "G4VPhysicalVolume.hh"


// Shared singleton -- instantiated during geom building on master only

G4CMPLogicalBorderTable* G4CMPLogicalBorderSurface::surfaceTable = 0;

const G4CMPLogicalBorderTable* G4CMPLogicalBorderSurface::GetSurfaceTable() {
  if (!surfaceTable) surfaceTable = new G4CMPLogicalBorderTable;
  return surfaceTable;
}


// Constructor

G4CMPLogicalBorderSurface::
G4CMPLogicalBorderSurface(const G4String& name,
			  G4VPhysicalVolume* vol1, G4VPhysicalVolume* vol2,
			  G4SurfaceProperty* surfaceProperty)
  : G4LogicalSurface(name, surfaceProperty), Volume1(0), Volume2(0) {
  if (!surfaceTable) surfaceTable = new G4CMPLogicalBorderTable;
  SetVolumes(vol1, vol2);		// Also registers into map
}


// Reassign object in lookup table

void G4CMPLogicalBorderSurface::SetVolumes(G4VPhysicalVolume* vol1,
					   G4VPhysicalVolume* vol2) {
  if (!vol1 || !vol2) {
    G4Exception("G4CMPLogicalBorderSurface", "G4CMPSurf002",
		FatalException, "Null pointer passed for placement volume.");
    return;
  }

  RemoveFromTable();
  Volume1 = vol1;
  Volume2 = vol2;
  AddToTable();
}

// Used by SetVolumes() for entering or replacing object in registry

void G4CMPLogicalBorderSurface::RemoveFromTable() {
  if (surfaceTable)
    surfaceTable->erase(G4CMPLogicalBorderKey(Volume1,Volume2));
}

void G4CMPLogicalBorderSurface::AddToTable() {
  if (!surfaceTable) GetSurfaceTable();
  (*surfaceTable)[G4CMPLogicalBorderKey(Volume1,Volume2)] = this;
}


// Access lookup table mapping surfaces to volumes

G4CMPLogicalBorderSurface* 
G4CMPLogicalBorderSurface::GetSurface(const G4VPhysicalVolume* vol1,
				      const G4VPhysicalVolume* vol2) {
  if (!surfaceTable) return 0;

  G4CMPLogicalBorderKey key(vol1,vol2);
  G4CMPLogicalBorderTable::iterator entry = surfaceTable->find(key);
  return (entry != surfaceTable->end() ? entry->second : 0);
}

void G4CMPLogicalBorderSurface::CleanSurfaceTable() {
  if (!surfaceTable) return;		// No table defined, nothing to clear

  // Delete all defined surfaces, then remove table entries
  for (auto pos=surfaceTable->cbegin(); pos!=surfaceTable->cend(); ++pos) {
    if (pos->second) delete pos->second;
  }

  surfaceTable->clear();
}


// Report contents of registry for diagnostics

void G4CMPLogicalBorderSurface::DumpInfo() {
  G4cout << "***** G4CMPLogicalBorderSurface Table : Nb of Surfaces = "
         << GetNumberOfSurfaces() << " *****" << G4endl;

  if (!surfaceTable) return;

  for (auto pos=surfaceTable->cbegin(); pos!=surfaceTable->cend(); ++pos) {
    G4cout << pos->second->GetName() << " : " << G4endl
	   << " Border of volumes "
	   << pos->second->GetVolume1()->GetName() << " and " 
	   << pos->second->GetVolume2()->GetName()
	   << G4endl;
  }

  G4cout << G4endl;
}

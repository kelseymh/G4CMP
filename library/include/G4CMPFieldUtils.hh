/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef G4CMPFieldUtils_hh
#define G4CMPFieldUtils_hh 1

// $Id$
// File: G4CMPFieldUtils.hh
//
// Description: Free standing helper functions for electric field access
//
// 20180622  Michael Kelsey
// 20211005  Add position-only utility (get touchable from position)

#include "G4ThreeVector.hh"

class G4CMPMeshElectricField;
class G4CMPLocalElectroMagField;
class G4LogicalVolume;
class G4Step;
class G4Track;
class G4UniformElectricField;
class G4VTouchable;


namespace G4CMP {
  // Get field at _global_ coordinate of track/step
  G4ThreeVector GetFieldAtPosition(const G4Track& track);
  G4ThreeVector GetFieldAtPosition(const G4Step& step);
  G4ThreeVector GetFieldAtPosition(const G4ThreeVector& pos);
  G4ThreeVector GetFieldAtPosition(const G4VTouchable* touch,
				   const G4ThreeVector& pos);
  G4ThreeVector GetFieldAtPosition(const G4LogicalVolume* vol,
				   const G4ThreeVector& pos);

  // Get potential (uniform or mesh field) at _global_ coordinate
  G4double GetPotentialAtPosition(const G4Step& step);
  G4double GetPotentialAtPosition(const G4Track& track);
  G4double GetPotentialAtPosition(const G4ThreeVector& pos);
  G4double GetPotentialAtPosition(const G4VTouchable* touch,
				  const G4ThreeVector& pos);
  G4double GetPotentialAtPosition(const G4LogicalVolume* vol,
				  const G4ThreeVector& pos);

  G4double GetPotentialAtVertex(const G4Track& track);

  // Estimate bias across volume through _global_ coordinate
  // NOTE:  Doesn't integrate field, just uses local magnitude, direction
  G4double GetBiasThroughPosition(const G4VTouchable* touch,
				  const G4ThreeVector& pos);

  // Return field handler if available for volume
  const G4CMPLocalElectroMagField* GetLocalField(const G4LogicalVolume* vol);
  const G4CMPMeshElectricField* GetMeshField(const G4LogicalVolume* vol);
  const G4UniformElectricField* GetUniformField(const G4LogicalVolume* vol);
}

#endif	/* G4CMPFieldUtils_hh */

/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
//
/// \file  library/src/G4CMPVElectrodePattern.cc
/// \brief Abstract base class to define complex electrode layouts
//

#include "G4CMPVElectrodePattern.hh"
#include "G4MaterialPropertiesTable.hh"


// Handles casting table to non-const for access

G4double
G4CMPVElectrodePattern::GetMaterialProperty(const G4String& key) const {
  G4MaterialPropertiesTable* ncTable =
    const_cast<G4MaterialPropertiesTable*>(theSurfaceTable);

  return ncTable->GetConstProperty(key);
}

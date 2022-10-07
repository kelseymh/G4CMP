/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
//
/// \file  library/src/G4CMPVElectrodePattern.cc
/// \brief Abstract base class to define complex electrode layouts
//
// 20221006  G4CMP-330 -- Pass lattice temperature through to sensors

#include "G4CMPVElectrodePattern.hh"
#include "G4LatticePhysical.hh"
#include "G4MaterialPropertiesTable.hh"


// Transfer lattice temperature into properties table if not already set

void G4CMPVElectrodePattern::
UseSurfaceTable(G4MaterialPropertiesTable* surfProp) {
  theSurfaceTable = surfProp;

  if (theLattice && !theSurfaceTable->ConstPropertyExists("temperature"))
    theSurfaceTable->AddConstProperty("temperature",
				      theLattice->GetTemperature());
}


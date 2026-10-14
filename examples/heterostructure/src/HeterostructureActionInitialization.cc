/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// /file  HeterostructureActionInitialization.cc
/// /brief Implementation of the HeterostructureActionInitialization class

// 20260815  Selby Q. Dang (Stanford/SLAC)

#include "HeterostructureActionInitialization.hh"
#include "HeterostructurePrimaryGeneratorAction.hh"
#include "G4CMPStackingAction.hh"

void HeterostructureActionInitialization::Build() const {
  SetUserAction(new HeterostructurePrimaryGeneratorAction);
  SetUserAction(new G4CMPStackingAction);
} 

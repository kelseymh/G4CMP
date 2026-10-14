/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file PhononActionInitialization.cc
/// \brief Implementation of the PhononActionInitialization class.
///
/// Creates PrimaryGenerator and StackingAction instances, following Geant4
/// requirements.  Additional UserAction classes should be added as needed.

#include "PhononActionInitialization.hh"
#include "PhononPrimaryGeneratorAction.hh"
#include "G4CMPStackingAction.hh"

void PhononActionInitialization::Build() const {
  SetUserAction(new PhononPrimaryGeneratorAction);
  SetUserAction(new G4CMPStackingAction);
} 

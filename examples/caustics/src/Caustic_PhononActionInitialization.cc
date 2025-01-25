/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// 20241024 Israel Hernandez -- IIT, QSC and Fermilab

#include "Caustic_PhononActionInitialization.hh"
#include "Caustic_PhononPrimaryGeneratorAction.hh"
#include "G4CMPStackingAction.hh"

void Caustic_PhononActionInitialization::Build() const {
  SetUserAction(new Caustic_PhononPrimaryGeneratorAction);
  SetUserAction(new G4CMPStackingAction);
}

/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#include "ChargeActionInitialization.hh"
#include "ChargePrimaryGeneratorAction.hh"
#include "G4CMPStackingAction.hh"

void ChargeActionInitialization::Build() const {
  SetUserAction(new ChargePrimaryGeneratorAction);
  SetUserAction(new G4CMPStackingAction);
} 

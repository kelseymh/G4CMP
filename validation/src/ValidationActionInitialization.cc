/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id: 539f524339ae53ad098a07cfa3bebd07784d23dd $

#include "ValidationActionInitialization.hh"
#include "ValidationPrimaryGeneratorAction.hh"
#include "ValidationSteppingAction.hh"
#include "G4CMPStackingAction.hh"


void ValidationActionInitialization::Build() const {
  SetUserAction(new ValidationPrimaryGeneratorAction);
  SetUserAction(new G4CMPStackingAction);
  SetUserAction(new ValidationSteppingAction);
} 

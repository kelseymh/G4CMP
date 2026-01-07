/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id: 539f524339ae53ad098a07cfa3bebd07784d23dd $

#include "QuasiparticleActionInitialization.hh"
#include "QuasiparticlePrimaryGeneratorAction.hh"
#include "QuasiparticleSteppingAction.hh"
#include "G4CMPStackingAction.hh"

void QuasiparticleActionInitialization::Build() const {
  SetUserAction(new QuasiparticlePrimaryGeneratorAction);
  SetUserAction(new G4CMPStackingAction);
  SetUserAction(new QuasiparticleSteppingAction);
} 

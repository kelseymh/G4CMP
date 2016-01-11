// $Id: 539f524339ae53ad098a07cfa3bebd07784d23dd $

#include "PhononActionInitialization.hh"
#include "PhononPrimaryGeneratorAction.hh"
#include "G4CMPStackingAction.hh"

void PhononActionInitialization::Build() const {
  SetUserAction(new PhononPrimaryGeneratorAction);
  SetUserAction(new G4CMPStackingAction);
} 

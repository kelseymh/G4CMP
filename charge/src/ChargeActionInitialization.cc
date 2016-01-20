#include "ChargeActionInitialization.hh"
#include "ChargePrimaryGeneratorAction.hh"
#include "G4CMPStackingAction.hh"

void ChargeActionInitialization::Build() const {
  SetUserAction(new ChargePrimaryGeneratorAction);
  SetUserAction(new G4CMPStackingAction);
} 

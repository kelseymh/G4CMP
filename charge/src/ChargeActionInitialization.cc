// $Id$
//
// For Geant4 10.0, replaces runManager->SetUserAction() calls.
// Not valid for previous releases; G4VERSION check is used to disable.

#include "ChargeActionInitialization.hh"
#include "ChargePrimaryGeneratorAction.hh"
#include "G4CMPStackingAction.hh"

#include "G4Version.hh"
#if (G4VERSION_NUMBER >= 1000)

void ChargeActionInitialization::Build() const {
  SetUserAction(new ChargePrimaryGeneratorAction);
  SetUserAction(new G4CMPStackingAction);
} 

#endif	/* G4VERSION_NUMBER */

// $Id$
//
// For Geant4 10.0, replaces runManager->SetUserAction() calls.
// Not valid for previous releases; G4VERSION check is used to disable.

#ifndef PhononActionInitialization_hh
#define PhononActionInitialization_hh 1

#include "G4Version.hh"
#if (G4VERSION_NUMBER >= 1000)

#include "G4VUserActionInitialization.hh"

class PhononActionInitialization : public G4VUserActionInitialization {
public:
  PhononActionInitialization() {;}
  virtual ~PhononActionInitialization() {;}
  virtual void Build() const;
};

#endif	/* G4VERSION_NUMBER */
#endif	/* PhononActionInitialization_hh */

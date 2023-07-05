/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// 20170525  Remove empty default destructor ("rule of five" semantics)

#ifndef G4CMPDriftElectron_h
#define G4CMPDriftElectron_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"

class G4CMPDriftElectron : public G4ParticleDefinition {
public:
  static G4CMPDriftElectron* Definition();
  static G4CMPDriftElectron* G4CMPDriftElectronDefinition();

private:
  static G4CMPDriftElectron* theInstance;
  G4CMPDriftElectron() {;}
};

#endif	/* G4CMPDriftElectron_h */

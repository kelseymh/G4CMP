/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// 20170525  Remove empty default destructor ("rule of five" semantics)

#ifndef G4CMPDriftHole_h
#define G4CMPDriftHole_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"

class G4CMPDriftHole : public G4ParticleDefinition {
public:
  static G4CMPDriftHole* Definition();
  static G4CMPDriftHole* G4CMPDriftHoleDefinition();

private:
  static G4CMPDriftHole* theInstance;
  G4CMPDriftHole() {;}
};

#endif	/* G4CMPDriftHole_h */

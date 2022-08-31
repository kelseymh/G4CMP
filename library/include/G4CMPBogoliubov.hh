/***********************************************************************\
 *  * This software is licensed under the terms of the GNU General Public *
 *   * License version 3 or later. See G4CMP/LICENSE for the full license. *
 *   \***********************************************************************/

///
///// \brief Definition of the G4CMPBogoliubovon class for quasiparticle tracking
////
//// $Id$
////

#ifndef G4CMPBogoliubov_h
#define G4CMPBogoliubov_h 1

#include "globals.hh"
#include "G4ParticleDefinition.hh"

class G4CMPBogoliubov : public G4ParticleDefinition {
public:
  static G4CMPBogoliubov* Definition();
  static G4CMPBogoliubov* G4CMPBogoliubovDefinition();

private:
  static G4CMPBogoliubov* theInstance;
  G4CMPBogoliubov() {;}
};

#endif	/* G4CMPBogoliubovon */

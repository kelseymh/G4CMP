/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id: e58a61fedbb99b167e16dafebc9c8664ae0c7b94 $

#ifndef QuasiparticleActionInitialization_hh
#define QuasiparticleActionInitialization_hh 1

#include "G4VUserActionInitialization.hh"

class QuasiparticleActionInitialization : public G4VUserActionInitialization {
public:
  QuasiparticleActionInitialization() {;}
  virtual ~QuasiparticleActionInitialization() {;}
  virtual void Build() const;
};

#endif	/* QuasiparticleActionInitialization_hh */

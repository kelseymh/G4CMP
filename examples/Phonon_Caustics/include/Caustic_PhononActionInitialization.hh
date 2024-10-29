/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// 20241024 Israel Hernandez -- IIT, QSC and Fermilab

#ifndef Caustic_PhononActionInitialization_hh
#define Caustic_PhononActionInitialization_hh 1

#include "G4VUserActionInitialization.hh"

class Caustic_PhononActionInitialization : public G4VUserActionInitialization {
public:
  Caustic_PhononActionInitialization() {;}
  virtual ~Caustic_PhononActionInitialization() {;}
  virtual void Build() const;
};

#endif

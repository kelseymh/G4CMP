/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file PhononActionInitialization.hh
/// \brief Definition of the PhononActionInitialization class.

#ifndef PhononActionInitialization_hh
#define PhononActionInitialization_hh 1

#include "G4VUserActionInitialization.hh"

class PhononActionInitialization : public G4VUserActionInitialization {
public:
  PhononActionInitialization() {;}
  virtual ~PhononActionInitialization() {;}
  virtual void Build() const;
};

#endif	/* PhononActionInitialization_hh */

/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file ValidationActionInitialization.hh
/// \brief Definition of the class governing action initialization

#ifndef ValidationActionInitialization_hh
#define ValidationActionInitialization_hh 1

#include "G4VUserActionInitialization.hh"

class ValidationActionInitialization : public G4VUserActionInitialization {
public:
  ValidationActionInitialization() {;}
  virtual ~ValidationActionInitialization() {;}
  virtual void Build() const;
};

#endif	/* ValidationActionInitialization_hh */

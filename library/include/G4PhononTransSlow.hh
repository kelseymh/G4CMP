/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file particles/phonons/include/G4PhononTransSlow.hh
/// \brief Definition of the G4PhononTransSlow class
//
// $Id$
//

#ifndef G4PhononTransSlow_h
#define G4PhononTransSlow_h 1

#include "G4ParticleDefinition.hh"

class G4PhononTransSlow : public G4ParticleDefinition {
private:
  static G4PhononTransSlow* theInstance;

private:
  G4PhononTransSlow() {;}

public:
  virtual ~G4PhononTransSlow () {;}
  
  static G4PhononTransSlow* Definition();
  static G4PhononTransSlow* PhononDefinition();
};

#endif	/* G4PhononTransSlow_h */

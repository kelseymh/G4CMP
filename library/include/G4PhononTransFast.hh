/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file particles/phonons/include/G4PhononTransFast.hh
/// \brief Definition of the G4PhononTransFast class
//
// $Id$
//

#ifndef G4PhononTransFast_h
#define G4PhononTransFast_h 1

#include "G4ParticleDefinition.hh"


class G4PhononTransFast : public G4ParticleDefinition {
private:
  static G4PhononTransFast* theInstance;
  
private:
  G4PhononTransFast() {;}
  
public:
  virtual ~G4PhononTransFast() {;}
  
  static G4PhononTransFast* Definition();
  static G4PhononTransFast* PhononDefinition();
};

#endif	/* G4PhononTransFast_h */

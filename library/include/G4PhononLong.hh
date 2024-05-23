/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file particles/include/G4PhononLong.hh
/// \brief Definition of the G4PhononLong class
//
// $Id$
//
#ifndef G4PhononLong_h
#define G4PhononLong_h 1

#include "G4ParticleDefinition.hh"


class G4PhononLong : public G4ParticleDefinition {
private:
  static G4PhononLong* theInstance;

private:
  G4PhononLong() {;}

public:
  virtual ~G4PhononLong() {;}

  static G4PhononLong* Definition();
  static G4PhononLong* PhononDefinition();

};

#endif	/* G4PhononLong_h */

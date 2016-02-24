/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4PhononPolarization.hh
/// \brief Definition of the G4PhononPolarization enum
//
// $Id$
//
#ifndef G4PhononPolarization_h
#define G4PhononPolarization_h 1

#include "globals.hh"

class G4ParticleDefinition;


namespace G4PhononPolarization {
  enum { Long=0, TransSlow=1, TransFast=2, UNKNOWN=-1 };

  G4int Get(const G4ParticleDefinition* aPD);
  G4ParticleDefinition* Get(G4int pol);
}

#endif	/* G4PhononPolarization_h */

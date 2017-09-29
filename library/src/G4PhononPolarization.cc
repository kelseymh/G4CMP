/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4PhononPolarization.cc
/// \brief implementation of the G4PhononPolarization enum
//
// $Id$
//
// 20160622  Add functions to return name strings

#include "G4PhononPolarization.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhononLong.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononTransSlow.hh"


G4int G4PhononPolarization::Get(const G4ParticleDefinition* aPD) {
  if (aPD == G4PhononLong::Definition())      return Long;
  if (aPD == G4PhononTransSlow::Definition()) return TransSlow;
  if (aPD == G4PhononTransFast::Definition()) return TransFast;
  return UNKNOWN;
}

G4ParticleDefinition* G4PhononPolarization::Get(G4int mode) {
  switch (mode) {
  case Long:      return G4PhononLong::Definition(); break;
  case TransSlow: return G4PhononTransSlow::Definition(); break;
  case TransFast: return G4PhononTransFast::Definition(); break;
  default: ;
  }

  return 0;
}

// Returns full particle name (convenience funciton)
const G4String& G4PhononPolarization::Name(G4int mode) {
  static const G4String empty;
  return (mode>=Long && mode<NUM_MODES) ? Get(mode)->GetParticleName() : empty;
}

// Returns short affix (L, ST, FT)
const char* G4PhononPolarization::Label(G4int mode) {
  switch (mode) {
  case Long:      return "L"; break;
  case TransSlow: return "ST"; break;
  case TransFast: return "LT"; break;
  default: ;
  }

  return "";
}

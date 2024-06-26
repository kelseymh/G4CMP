/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file particles/include/G4BogoliubovQP.hh
/// \brief Definition of the G4BogoliubovQP class
// This class deals with quasiparticles in the superconducting layer just as the Kaplan
// class did, but now with more rigorous Geant4 particle tracking. 


#ifndef G4BogoliubovQP_h
#define G4BogoliubovQP_h 1

#include "G4ParticleDefinition.hh"


class G4BogoliubovQP : public G4ParticleDefinition {
private:
  static G4BogoliubovQP* theInstance;

private:
  G4BogoliubovQP() {;}

public:
  virtual ~G4BogoliubovQP() {;}

  static G4BogoliubovQP* Definition();
  static G4BogoliubovQP* BogoliubovQPDefinition();

};

#endif	/* G4BogoliubovQP_h */


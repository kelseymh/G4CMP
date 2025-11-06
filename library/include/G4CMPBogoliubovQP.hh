/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file particles/include/G4CMPBogoliubovQP.hh
/// \brief Definition of the G4CMPBogoliubovQP class
///
/// This class deals with quasiparticles in the superconducting layer
/// just as the Kaplan class did, but now with more rigorous Geant4
/// particle tracking.
//
// 20250922  G4CMP-219 -- First addition to this history (done at time
//                        of merge to develop)

#ifndef G4CMPBogoliubovQP_h
#define G4CMPBogoliubovQP_h 1

#include "G4ParticleDefinition.hh"


class G4CMPBogoliubovQP : public G4ParticleDefinition {
private:
  static G4CMPBogoliubovQP* theInstance;

private:
  G4CMPBogoliubovQP() {;}

public:
  virtual ~G4CMPBogoliubovQP() {;}

  static G4CMPBogoliubovQP* Definition();
  static G4CMPBogoliubovQP* BogoliubovQPDefinition();
};

#endif	/* G4CMPBogoliubovQP_h */

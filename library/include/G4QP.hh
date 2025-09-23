/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file particles/include/G4QP.hh
/// \brief Definition of the G4QP class
///
/// This class deals with quasiparticles in the superconducting layer
/// just as the Kaplan class did, but now with more rigorous Geant4
/// particle tracking.
//
// 20250922  G4CMP-219 -- First addition to this history (done at time
//                        of merge to develop)

#ifndef G4QP_h
#define G4QP_h 1

#include "G4ParticleDefinition.hh"


class G4QP : public G4ParticleDefinition {
private:
  static G4QP* theInstance;

private:
  G4QP() {;}

public:
  virtual ~G4QP() {;}

  static G4QP* Definition();
  static G4QP* QPDefinition();
};

#endif	/* G4QP_h */

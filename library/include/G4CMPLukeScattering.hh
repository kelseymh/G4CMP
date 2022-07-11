/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPLukeScattering.hh
/// \brief Definition of the G4CMPLukeScattering class
//
// $Id$
//
// 20150111  New base class for both electron and hole Luke processes
// 20160110  Remerge the electron and hole subclasses into one class
// 20170805  Remove GetMeanFreePath() function to scattering-rate model
// 20190816  Add flag to track secondary phonons immediately (c.f. G4Cerenkov)
// 20201109  Drop G4CMP_DEBUG protection here, to avoid client rebuilding

#ifndef G4CMPLukeScattering_h
#define G4CMPLukeScattering_h 1

#include "globals.hh"
#include "G4CMPVDriftProcess.hh"
#include "G4ThreeVector.hh"
#include <iostream>

class G4CMPTrackInformation;
class G4VProcess;
class G4ParticleDefinition;
class G4Track;

class G4CMPLukeScattering : public G4CMPVDriftProcess {
public:
  G4CMPLukeScattering(G4VProcess* stepper=0);
  virtual ~G4CMPLukeScattering();

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

  // Pause current particle tracking, track secondary phonons instead
  void SetTrackSecondariesFirst(const G4bool val) { secondariesFirst = val; }
  G4bool GetTrackSecondariesFirst() const { return secondariesFirst; }

private:
  // hide assignment operator as private
  G4CMPLukeScattering(G4CMPLukeScattering&);
  G4CMPLukeScattering& operator=(const G4CMPLukeScattering& right);

  G4VProcess* stepLimiter;
  G4bool secondariesFirst;

  std::ofstream output;		// Only used for G4CMP_DEBUG debugging
};

#endif	/* G4CMPLukeScattering */

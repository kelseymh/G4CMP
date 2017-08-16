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

private:
  // hide assignment operator as private
  G4CMPLukeScattering(G4CMPLukeScattering&);
  G4CMPLukeScattering& operator=(const G4CMPLukeScattering& right);

  G4VProcess* stepLimiter;
#ifdef G4CMP_DEBUG
  std::ofstream output;
#endif
};

#endif	/* G4CMPLukeScattering */

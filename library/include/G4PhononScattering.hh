/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4PhononScattering.hh
/// \brief Definition of the G4PhononScattering class
//
// $Id$
//
// 20170805  Move GetMeanFreePath() to scattering-rate model

#ifndef G4PhononScattering_h
#define G4PhononScattering_h 1

#include "G4VPhononProcess.hh"

class G4PhononScattering : public G4VPhononProcess {
public:
  G4PhononScattering(const G4String& processName="phononScattering");
  virtual ~G4PhononScattering();

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

protected:
  // Keep function here as call-back to avoid getting old toolkit version
  virtual G4double GetMeanFreePath(const G4Track& trk, G4double prevstep,
				   G4ForceCondition* cond) {
    return G4CMPVProcess::GetMeanFreePath(trk, prevstep, cond);
  }

private:
  // hide assignment operator as private 
  G4PhononScattering(G4PhononScattering&);
  G4PhononScattering& operator=(const G4PhononScattering& right);
};

#endif	/* G4PhononScattering_h */

/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPQPLocalTrappingProcess.hh
/// \brief Definition of the G4CMPQPLocalTrappingProcess class
//
//
// 20170805  Move GetMeanFreePath() to scattering-rate model

#ifndef G4CMPQPLocalTrappingProcess_h
#define G4CMPQPLocalTrappingProcess_h 1

#include "G4VQPProcess.hh"

class G4CMPQPLocalTrappingProcess : public G4VQPProcess {
public:
  G4CMPQPLocalTrappingProcess(const G4String& processName="QPLocalTrapping");
  virtual ~G4CMPQPLocalTrappingProcess();

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

protected:
  // Keep function here as call-back to avoid getting old toolkit version
  virtual G4double GetMeanFreePath(const G4Track& trk,
				   G4double prevstep,
				   G4ForceCondition* cond);

private:
  // hide assignment operator as private 
  G4CMPQPLocalTrappingProcess(G4CMPQPLocalTrappingProcess&);
  G4CMPQPLocalTrappingProcess& operator=(const G4CMPQPLocalTrappingProcess& right);
};

#endif	/* G4CMPQPLocalTrappingProcess_h */

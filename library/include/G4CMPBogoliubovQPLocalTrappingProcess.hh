/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPBogoliubovQPLocalTrappingProcess.hh
/// \brief Definition of the G4CMPBogoliubovQPLocalTrappingProcess class
//
//
// 20170805  Move GetMeanFreePath() to scattering-rate model

#ifndef G4CMPBogoliubovQPLocalTrappingProcess_h
#define G4CMPBogoliubovQPLocalTrappingProcess_h 1

#include "G4VBogoliubovQPProcess.hh"

class G4CMPBogoliubovQPLocalTrappingProcess : public G4VBogoliubovQPProcess {
public:
  G4CMPBogoliubovQPLocalTrappingProcess(const G4String& processName="bogoliubovQPLocalTrapping");
  virtual ~G4CMPBogoliubovQPLocalTrappingProcess();

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

protected:
  // Keep function here as call-back to avoid getting old toolkit version
  virtual G4double GetMeanFreePath(const G4Track& trk,
				   G4double prevstep,
				   G4ForceCondition* cond);

private:
  // hide assignment operator as private 
  G4CMPBogoliubovQPLocalTrappingProcess(G4CMPBogoliubovQPLocalTrappingProcess&);
  G4CMPBogoliubovQPLocalTrappingProcess& operator=(const G4CMPBogoliubovQPLocalTrappingProcess& right);
};

#endif	/* G4CMPBogoliubovQPLocalTrappingProcess_h */

/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPQPLocalTrappingProcess.hh
/// \brief Definition of the G4CMPQPLocalTrappingProcess class
///
/// This process captures QPs' ability to trap on local trapping sites
/// due to (non-engineered) local fluctuations in the superconducting
/// gap or impurities. It is described by a simple trapping time.
//
// 20250922  G4CMP-219 -- First addition to this history (done at time
//                        of merge to develop)

#ifndef G4CMPQPLocalTrappingProcess_h
#define G4CMPQPLocalTrappingProcess_h 1

#include "G4CMPVQPProcess.hh"

class G4CMPQPLocalTrappingProcess : public G4CMPVQPProcess {
public:
  G4CMPQPLocalTrappingProcess(const G4String& processName="qpLocalTrapping");
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

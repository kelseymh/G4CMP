/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPQPDiffusionTimeStepperProcess.hh
/// \brief Definition of the G4CMPQPDiffusionTimeStepperProcess class
//

#ifndef G4CMPQPDiffusionTimeStepperProcess_h
#define G4CMPQPDiffusionTimeStepperProcess_h 1

#include "G4CMPSCUtils.hh"
#include "G4VQPProcess.hh"

class G4CMPQPDiffusionTimeStepperProcess : public G4VQPProcess
{
public:
  G4CMPQPDiffusionTimeStepperProcess(const G4String& processName ="qpDiffusionTimeStepper");
  virtual ~G4CMPQPDiffusionTimeStepperProcess();

  virtual void SetVerboseLevel(G4int vb);
  virtual G4bool IsApplicable(const G4ParticleDefinition&);
  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step& );

protected:
  // Keep function here as call-back to avoid getting old toolkit version
  virtual G4double GetMeanFreePath(const G4Track& trk,
				   G4double prevstep,
				   G4ForceCondition* cond);



private:

  // hide assignment operator as private
  G4CMPQPDiffusionTimeStepperProcess(G4CMPQPDiffusionTimeStepperProcess&);
  G4CMPQPDiffusionTimeStepperProcess& operator=(const G4CMPQPDiffusionTimeStepperProcess& right);
};

#endif	/* G4CMPQPDiffusionTimeStepperProcess_h */

/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPQPDiffusionTimeStepperProcess.hh
/// \brief Definition of the G4CMPQPDiffusionTimeStepperProcess class
///
/// This process is an auxiliary helper process for QP diffusion. It is
/// primarily intended as a debugging tool. It is a "dummy" process
/// that can be used to compete with the main QP diffusion process
/// which uses the Walk-on-Spheres approach. This one just simply
/// assumes a characteristic timescale for the next diffusion "scatter"
/// that, unlike the WoS approach, is agnostic to the 2D safety.
/// So this acts a bit like all of the other inelastic QP processes
/// as they relate to the main QP diffusion process.
///
/// A few extra notes/helpful tidbits:
/// 1. QP diffusion should run perfectly well if this is process is
///    turned off. It should be turned off by default.
/// 2. If this process is turned on, diffusion should in principle
///    yield the same underlying physics, but the set of steps that
///    a QP takes may be different depending on the value of the
///    time stepper length used. Higher time stepper rates will
///    increase the number and granularity of steps between "real"
///    (i.e. inelastic) processes undergone by the QP.
/// 3. The parameter for this time stepper process is hardcoded into
///    the time stepper rate .cc file. It is deliberately done so
///    so that people don't feel encouraged to muck with that rate
///    given that this is primarily intended as a debugging tool and
///    given that the main QP diffusion function should be perfectly
///    fine without this process helping.
//
// 20250922  G4CMP-219 : First addition to this history (done at time of
//           merge to develop)


#ifndef G4CMPQPDiffusionTimeStepperProcess_h
#define G4CMPQPDiffusionTimeStepperProcess_h 1

#include "G4CMPSCUtils.hh"
#include "G4VQPProcess.hh"

class G4CMPQPDiffusionTimeStepperProcess : public G4VQPProcess {
public:
  G4CMPQPDiffusionTimeStepperProcess(const G4String& processName ="qpDiffusionTimeStepper");
  virtual ~G4CMPQPDiffusionTimeStepperProcess();

  virtual void SetVerboseLevel(G4int vb);
  virtual G4bool IsApplicable(const G4ParticleDefinition&);
  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step& );

protected:
  // Keep function here as call-back to avoid getting old toolkit version
  virtual G4double GetMeanFreePath(const G4Track& trk,G4double prevstep,
				   G4ForceCondition* cond);



private:

  // hide assignment operator as private
  G4CMPQPDiffusionTimeStepperProcess(G4CMPQPDiffusionTimeStepperProcess&);
  G4CMPQPDiffusionTimeStepperProcess& operator=(const G4CMPQPDiffusionTimeStepperProcess& right);
};

#endif	/* G4CMPQPDiffusionTimeStepperProcess_h */

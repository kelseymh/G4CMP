/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPQPRecombinationProcess.hh
/// \brief Definition of the G4CMPQPRecombinationProcess class
//

#ifndef G4CMPQPRecombinationProcess_h
#define G4CMPQPRecombinationProcess_h 1

#include "G4CMPSCUtils.hh"
#include "G4VQPProcess.hh"

class G4CMPQPRecombinationProcess : public G4VQPProcess
{
public:
  G4CMPQPRecombinationProcess(const G4String& processName ="QPRecombination");
  virtual ~G4CMPQPRecombinationProcess();

  virtual void SetVerboseLevel(G4int vb);
  virtual G4bool IsApplicable(const G4ParticleDefinition&);
  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step& );

protected:
  // Keep function here as call-back to avoid getting old toolkit version
  virtual G4double GetMeanFreePath(const G4Track& trk, G4double prevstep,
				   G4ForceCondition* cond);



private:

  void GenerateRecombinationPhonon(G4double phonEnergy,
				   const G4Track& aTrack,
				   const G4Step& aStep);
  

  // hide assignment operator as private
  G4CMPQPRecombinationProcess(G4CMPQPRecombinationProcess&);
  G4CMPQPRecombinationProcess& operator=(const G4CMPQPRecombinationProcess& right);
};

#endif	/* G4CMPQPRecombinationProcess_h */

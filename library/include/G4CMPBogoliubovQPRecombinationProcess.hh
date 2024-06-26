/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPBogoliubovQPRecombinationProcess.hh
/// \brief Definition of the G4CMPBogoliubovQPRecombinationProcess class
//

#ifndef G4CMPBogoliubovQPRecombinationProcess_h
#define G4CMPBogoliubovQPRecombinationProcess_h 1

#include "G4CMPSCUtils.hh"
#include "G4VBogoliubovQPProcess.hh"

class G4CMPBogoliubovQPRecombinationProcess : public G4VBogoliubovQPProcess
{
public:
  G4CMPBogoliubovQPRecombinationProcess(const G4String& processName ="bogoliubovQPRecombination");
  virtual ~G4CMPBogoliubovQPRecombinationProcess();

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
  G4CMPBogoliubovQPRecombinationProcess(G4CMPBogoliubovQPRecombinationProcess&);
  G4CMPBogoliubovQPRecombinationProcess& operator=(const G4CMPBogoliubovQPRecombinationProcess& right);
};

#endif	/* G4CMPBogoliubovQPRecombinationProcess_h */

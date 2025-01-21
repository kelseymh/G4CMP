/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPBogoliubovQPRadiatesPhononProcess.hh
/// \brief Definition of the G4CMPBogoliubovQPRadiatesPhononProcess class
//

#ifndef G4CMPBogoliubovQPRadiatesPhononProcess_h
#define G4CMPBogoliubovQPRadiatesPhononProcess_h 1

#include "G4CMPSCUtils.hh"
#include "G4VBogoliubovQPProcess.hh"

class G4CMPBogoliubovQPRadiatesPhononProcess : public G4VBogoliubovQPProcess
{
public:
  G4CMPBogoliubovQPRadiatesPhononProcess(const G4String& processName ="bogoliubovQPRadiatesPhonon");
  virtual ~G4CMPBogoliubovQPRadiatesPhononProcess();

  virtual void SetVerboseLevel(G4int vb);
  virtual G4bool IsApplicable(const G4ParticleDefinition&);
  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step& );

protected:
  // Keep function here as call-back to avoid getting old toolkit version
  virtual G4double GetMeanFreePath(const G4Track& trk,
				   G4double prevstep,
				   G4ForceCondition* cond);



private:

  void GenerateRadiatedPhonon(G4double phonEnergy,
			      const G4Track& aTrack,
			      const G4Step& aStep);
  
  G4double PhononEnergyRand(G4double Energy) const;
  G4double PhononEnergyPDF(G4double E, G4double x) const;

  
  // hide assignment operator as private
  G4CMPBogoliubovQPRadiatesPhononProcess(G4CMPBogoliubovQPRadiatesPhononProcess&);
  G4CMPBogoliubovQPRadiatesPhononProcess& operator=(const G4CMPBogoliubovQPRadiatesPhononProcess& right);
};

#endif	/* G4CMPBogoliubovQPRadiatesPhononProcess_h */

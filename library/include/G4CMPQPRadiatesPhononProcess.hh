/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPQPRadiatesPhononProcess.hh
/// \brief Definition of the G4CMPQPRadiatesPhononProcess class
//

#ifndef G4CMPQPRadiatesPhononProcess_h
#define G4CMPQPRadiatesPhononProcess_h 1

#include "G4CMPSCUtils.hh"
#include "G4VQPProcess.hh"

class G4CMPQPRadiatesPhononProcess : public G4VQPProcess
{
public:
  G4CMPQPRadiatesPhononProcess(const G4String& processName ="qpRadiatesPhonon");
  virtual ~G4CMPQPRadiatesPhononProcess();

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
  G4CMPQPRadiatesPhononProcess(G4CMPQPRadiatesPhononProcess&);
  G4CMPQPRadiatesPhononProcess& operator=(const G4CMPQPRadiatesPhononProcess& right);
};

#endif	/* G4CMPQPRadiatesPhononProcess_h */

/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPSCPairBreakingProcess.hh
/// \brief Definition of the G4CMPSCPairBreakingProcess class
//

#ifndef G4CMPSCPairBreakingProcess_h
#define G4CMPSCPairBreakingProcess_h 1

#include "G4CMPSCUtils.hh"
#include "G4VPhononProcess.hh"

class G4CMPSCPairBreakingProcess : public G4VPhononProcess
{
public:
  G4CMPSCPairBreakingProcess(const G4String& processName ="scPairBreaking");
  virtual ~G4CMPSCPairBreakingProcess();

  virtual void SetVerboseLevel(G4int vb);
  virtual G4bool IsApplicable(const G4ParticleDefinition&);
  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step& );

protected:
  // Keep function here as call-back to avoid getting old toolkit version
  virtual G4double GetMeanFreePath(const G4Track& trk, G4double prevstep,
				   G4ForceCondition* cond);
  
  std::pair<G4double,G4double> FetchQPEnergies(G4double phonEnergy);
  G4double QPEnergyRand(G4double Energy) const;
  G4double QPEnergyPDF(G4double E, G4double x) const;

  
  void GenerateBogoliubovQPPair(std::pair<G4double,G4double> QPEnergies,
				const G4Track& aTrack,
				const G4Step& aStep);



private:

  

  // hide assignment operator as private
  G4CMPSCPairBreakingProcess(G4CMPSCPairBreakingProcess&);
  G4CMPSCPairBreakingProcess& operator=(const G4CMPSCPairBreakingProcess& right);
};

#endif	/* G4CMPSCPairBreakingProcess_h */

/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// 20200331  C. Stanford G4CMP-195: Added charge trapping
// 20200504  M. Kelsey: Reduce process name string; add static function to
//		return MFP by particle type

#ifndef G4CMPDriftTrappingProcess_h
#define G4CMPDriftTrappingProcess_h 1

#include "G4CMPVDriftProcess.hh"

class G4CMPEnergyPartition;


class G4CMPDriftTrappingProcess : public G4CMPVDriftProcess {
public:
  G4CMPDriftTrappingProcess(const G4String& name = "ChargeTrapping");
  virtual ~G4CMPDriftTrappingProcess();

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

  static G4double GetMeanFreePath(const G4ParticleDefinition* pd);

protected:
  virtual G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*);

private:
  G4CMPEnergyPartition* partitioner;

  // No copying/moving
  G4CMPDriftTrappingProcess(G4CMPDriftTrappingProcess&);
  G4CMPDriftTrappingProcess(G4CMPDriftTrappingProcess&&);
  G4CMPDriftTrappingProcess& operator=(const G4CMPDriftTrappingProcess&);
  G4CMPDriftTrappingProcess& operator=(const G4CMPDriftTrappingProcess&&);
};

#endif	/* G4CMPDriftTrappingProcess_h */

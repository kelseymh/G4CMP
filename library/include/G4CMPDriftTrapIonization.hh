/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// 20200331  G4CMP-196: Added impact ionization process
// 20200426  G4CMP-196: Change name to TrapIonization, specify beam and trap
//		particle types

#ifndef G4CMPDriftTrapIonization_h
#define G4CMPDriftTrapIonization_h 1

#include "G4CMPVDriftProcess.hh"

class G4ParticleDefinition;


class G4CMPDriftTrapIonization : public G4CMPVDriftProcess {
public:
  G4CMPDriftTrapIonization(G4ParticleDefinition* impactPD,
			   G4ParticleDefinition* trapPD,
			   const G4String& name="TrapIonization");
  virtual ~G4CMPDriftTrapIonization();

  // Separate process instances (separate MFPs) for each charge carrier
  virtual G4bool IsApplicable(const G4ParticleDefinition& aPD) {
    return (&aPD == impactType);
  }

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

  // Compute MFP for specific beam and trap types
  static G4double GetMeanFreePath(const G4ParticleDefinition* impactPD,
				  const G4ParticleDefinition* trapPD);

protected:
  virtual G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*);

protected:
  G4ParticleDefinition* impactType;	// e- or h+ for tracking
  G4ParticleDefinition* trapType;	// e- or h+ impurity trap

private:
  // No copying/moving
  G4CMPDriftTrapIonization(G4CMPDriftTrapIonization&);
  G4CMPDriftTrapIonization(G4CMPDriftTrapIonization&&);
  G4CMPDriftTrapIonization& operator=(const G4CMPDriftTrapIonization&);
  G4CMPDriftTrapIonization& operator=(const G4CMPDriftTrapIonization&&);
};

#endif	/* G4CMPDriftTrapIonization_h */
